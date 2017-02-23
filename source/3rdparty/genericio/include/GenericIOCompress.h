#ifndef GENERICIOCOMPRESS_H
#define GENERICIOCOMPRESS_H

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>
#include <string>
#include <cassert>

#include <vector>
#include <map>
#include <limits>
#include <algorithm>
#include <deque>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "AssocVector.h"

namespace gio {

// An implementation of Arithmetic Coding, see:
//
// Matt Mahoney. "Data Compression Explained." 2010-2013
// http://mattmahoney.net/dc/dce.html
class arithmetic_coding {
public:
  arithmetic_coding(std::vector<unsigned char> &d, bool encoding = true)
    : data(d), idx(0), low(0), high(0xFFFFFFFF), curr(0)
  {
    assert((encoding || data.size()) && "Decoding with no data?");

    if (!encoding) {
      for (; idx < 4; ++idx)
        curr = curr << 8 | (idx < data.size() ? data[idx] : 0);
    }
  }

  void encode(float q, int y) {
    assert((y == 0 || y == 1) && "Invalid bit value");
    assert(q >= 0 && q <= 1 && "Invalid probability value");
    assert(high > low && "Invalid state");

    // Represent probability q as an integer p as q = p/65535.
    int p = q*65535;
    uint32_t mid = low + (((uint64_t)(high - low) * p) >> 16);
    if (y)
      high = mid;
    else
      low = mid + 1;

    // Write out any identical leading bytes.
    while ((high^low) < 0x1000000) {
      data.push_back(high >> 24 /* same as low >> 24 */);
      high = (high << 8) | 255;
      low = low << 8;
    }
  }

  void flush() {
    assert(high > low && "Invalid state");
    data.push_back(high >> 24);
  }

  int decode(float q) {
    assert(q >= 0 && q <= 1 && "Invalid probability value");
    assert(high > low && "Invalid state");

    // Represent probability q as an integer p as q = p/65535.
    int p = q*65535;
    uint32_t mid = low + (((uint64_t)(high - low) * p) >> 16);

    int y = curr <= mid;
    if (y)
      high = mid;
    else
      low = mid + 1;

    // Shift out any identical leading bytes.
    while ((high^low) < 0x1000000) {
      high = (high << 8) | 255;
      low = low << 8;

      unsigned char c = 0;
      if (idx < data.size())
        c = data[idx++];
      curr = curr << 8 | c;
    }

    return y;
  }

protected:
  std::vector<unsigned char> &data;
  size_t idx; // current index for decoding.
  uint32_t low, high, curr;
};

template <bool B, typename T>
struct value_classification_ {
  static bool finite(T) { return true; }
  static bool nan(T) { return false; }
};

template <typename T>
struct value_classification_<true, T> {
  static bool finite(T v) { using namespace std; return isfinite(v); }
  static bool nan(T v) { using namespace std; return isnan(v); }
};

template <typename T>
struct value_classification
  : value_classification_<std::numeric_limits<T>::has_infinity  ||
                          std::numeric_limits<T>::has_quiet_NaN ||
                          std::numeric_limits<T>::has_signaling_NaN, T> {};

template <typename T>
class histogram_simple {
public:
  histogram_simple(size_t nb) : nBuckets(nb), ies(0.0) {}

  void add(T sample) {
    typename std::vector<T>::iterator i = lower_bound(elems.begin(), elems.end(), sample);
#if 0
    if (i == elems.end() || sample < *i)
      elems.insert(i, sample);
#endif
    elems.insert(i, sample);

    eque.push_back(sample);
    if (eque.size() > nBuckets) {
      T old = eque.front();
      i = find(elems.begin(), elems.end(), old);
      assert(i != elems.end() && "Cannot find old sample");
      elems.erase(i);
      eque.pop_front();
    }

#ifdef DEBUG_HISTOGRAM
    std::cout << "queue is now " << eque.size() << " elems is now " <<
            elems.size() << "\n";
#endif

    ies = 1.0f/float(elems.size());
  }

  // Estimate P(X <= x).
  float estimate_cdf(T x) const {
    if (!elems.size() || !value_classification<T>::finite(x))
      return 0.5f;

    typename std::vector<T>::const_iterator i = upper_bound(elems.begin(), elems.end(), x);
    size_t n = distance(elems.begin(), i);

#ifdef DEBUG_HISTOGRAM
    std::cout << "cdf of " << x << " " << n << " of " << elems.size() <<
            " elems not above\n";
#endif

    float q = float(n)*ies;
    if (i != elems.begin() && i != elems.end()) {
      q += ((x - *(i-1))/(*i - *(i-1)))/float(elems.size());
    }

    return q;
  }

  void print(std::ostream &os) const {}

protected:
  size_t nBuckets;
  std::vector<T> elems;
  std::deque<T> eque;

  // Dividing by float(elems.size()) is expensive, cache the inverse.
  float ies;
};

// A dynamic histogram based on:
//
// Donko Donjerkovic, Yannis Ioannidis, Raghu Ramakrishnan.
// "Dynamic Histograms: Capturing Evolving Data Sets."
// UW-Madison Technical Report CS-TR-99-1396, March, 1999
// http://pages.cs.wisc.edu/~donjerko/hist.pdf
template <typename T>
class histogram {
protected:
  struct bucket {
    bucket(size_t c = 0) : count(c), ucount(0), ccount(0), stamp(0) {}
    size_t count, ucount, ccount, stamp;
  };

  typedef AssocVector<T, bucket> map_type;
  typedef typename map_type::iterator map_iterator;
  typedef typename map_type::const_iterator const_map_iterator;

public:
  histogram(size_t nb) : nBuckets(nb), nTotal(0), maxAge(3000), curAge(0) {
    // To simplify the implementation, always start with a full-range bucket.
    // In ths way the total range of all numbers is covered.

    min_value = std::numeric_limits<T>::max();
    max_value = -std::numeric_limits<T>::max();

    buckets.insert(make_pair(std::numeric_limits<T>::has_infinity ?
                     std::numeric_limits<T>::infinity() : std::numeric_limits<T>::max(),
                   bucket()));
  }

  void add(T sample) {
    add_inner(sample);
    update_ccounts();
    ++curAge;
  }

  void print(std::ostream &os) const {
    size_t total = 0;
    for (const_map_iterator i = buckets.begin(), e = buckets.end();
         i != e; ++i) {
      T left = get_left_value(i);
      T right = get_right_value(i);
      os << left << "\t" << right << "\t" << i->second.count <<
#ifdef DEBUG_HISTOGRAM
        "\t" << i->second.ucount <<
#endif
        "\n";
      total += i->second.count;
    }

    os << "# buckets: " << buckets.size() << " total: " << total << "\n";
  }

  // Estimate P(X <= x).
  float estimate_cdf(T x) const {
    if (!nTotal || !value_classification<T>::finite(x))
      return 0.5f;

#ifdef DEBUG_HISTOGRAM
    std::cout << "compute cdf of x = " << x << "\n";
#endif

    const_map_iterator bi = buckets.lower_bound(x);
    assert(bi != buckets.end() && "Could not find bucket for sample");
    if (bi->first == x && bi != buckets.end()-1)
      ++bi;

    float counts = bi->second.ccount;

#ifdef DEBUG_HISTOGRAM
    std::cout << "found bucket " << bi->first << " for x = " << x << " with ccounts = " << counts << "\n";
#endif

    T left = get_left_value(bi);
    T right = get_right_value(bi);

#ifdef DEBUG_HISTOGRAM
    std::cout << "left = " << left << ", right = " << right << "\n";
#endif

    if (x >= left) {
      if (x < right) {
        // Capture the PDF contribution from x itself (which is taken here to
        // extend only to the right).
        T rx;
        if (std::numeric_limits<T>::is_integer)
          rx = x + T(1);
        else
          rx = x * (T(1) + std::numeric_limits<T>::epsilon());

        float range = (float) right - (float) left,
              part  = (float) rx    - (float) left;

#ifdef DEBUG_HISTOGRAM
        std::cout << "range = " << range << ", part = " << part << "\n";
#endif
        counts += range > 0 ? (part/range) * bi->second.count : 0;
      } else
        counts += bi->second.count;

#ifdef DEBUG_HISTOGRAM
      std::cout << "cdf counts: " << counts << "\n";
#endif
    }

#ifdef DEBUG_HISTOGRAM
    std::cout << "cdf of x = " << x << " is " << (counts / nTotal) << "\n";
#endif

    return counts / nTotal;
  }

protected:
  static bool bucket_is_zero(const std::pair<T, bucket> &p) {
    return !p.second.count;
  }

  void update_ccounts() {
    size_t c = 0;
    for (map_iterator i = buckets.begin(), e = buckets.end(); i != e; ++i) {
      i->second.ccount = c;
      c += i->second.count;
    }
  }

  void add_inner(T sample) {
    if (!value_classification<T>::finite(sample))
      return;
    map_iterator bi = (!std::numeric_limits<T>::has_infinity &&
                       sample == std::numeric_limits<T>::max()) ?
                      buckets.find(sample) : buckets.upper_bound(sample);
    assert(bi != buckets.end() && "Could not find bucket for sample");

#ifdef DEBUG_HISTOGRAM
    std::cout << "incrementing bucket " << bi->first << " from sample: " << sample << "\n";
#endif
    ++bi->second.count;
    ++nTotal;

    bi->second.stamp = curAge;

    min_value = std::min(min_value, sample);
    max_value = std::max(max_value, sample);

    T left = get_left_value(bi);
    if (!value_classification<T>::finite(left) ||
        sample > (bi->first/2 + left/2)) // distributed division to avoid overflows.
                                         // (note that if bi->first is infinite then
                                         //  the condition will be false).
      ++bi->second.ucount;

    if (buckets.size() < nBuckets) {
      if (bi->second.count == 1)
        return;

      // Don't split a zero-width bin.
      if (left == bi->first)
        return;

      // Don't split on the left value.
      if (sample == left)
        return;

#ifdef DEBUG_HISTOGRAM
      std::cout << "splitting bucket: " << bi->first << "\n";
#endif

      // We don't yet have the target number of buckets: split this one.
      split_bucket(bi, left, true /* use_sample */, sample);
      return;
    }

    float seps, meps;
    map_iterator s = best_to_split(seps, bi, sample), m = best_to_merge(meps),
                 e = buckets.end();
    if (seps > meps && s != e && s->second.count > 1 && m != e) {
#ifdef DEBUG_HISTOGRAM
      std::cout << "splitting bucket: " << s->first << " merging: " << m->first << "\n";
#endif

      T mfirst = m->first;
      bool use_sample = (s == buckets.begin() ||
                         s->first == buckets.rbegin()->first) && bi == s;
      split_bucket(s, get_left_value(s), use_sample, sample);

      // With a vector-based container, the split operation may invalidate the
      // m iterator. Find it again by value.
      m = buckets.find(mfirst);
      map_iterator ml = m;
      --ml;

      m->second.ucount = m->second.count;
      m->second.count += ml->second.count;
      m->second.stamp = max(m->second.stamp, ml->second.stamp);
      buckets.erase(ml);
    }

    for (map_iterator i = buckets.begin(), ie = buckets.end(); i != ie; ++i) {
      if (curAge - i->second.stamp > maxAge)
        i->second.count = i->second.ucount = 0;
    }

    // Remove zero bins (except the top-most bin, as it must always exist);
    map_iterator z = remove_if(buckets.begin(), buckets.end()-1, bucket_is_zero);
    e = buckets.end()-1;
    if (z != e)
      buckets.erase(z, e);
   
  }


  T get_left_value(const_map_iterator bi) const {
    if (bi == buckets.begin())
      return min_value;

    const_map_iterator bip = bi;
    --bip;
    return bip->first;
  }

  T get_right_value(const_map_iterator bi) const {
    if (bi->first == buckets.rbegin()->first) {
      if (std::numeric_limits<T>::is_integer)
        return max_value + T(1);
      else
        return max_value * (T(1) + std::numeric_limits<T>::epsilon());
    }

    return bi->first;
  }

  void split_bucket(map_iterator &bi, T left,
                    bool use_sample = false, T sample = T(0)) {
    bucket orig = bi->second;
    T right = bi->first;
    assert(left != right && "Splitting a single-value bin");

    // The original bucket becomes the new right-side bucket, and we insert a
    // new bucket to its left.
    bi->second.count = orig.ucount;

    bucket lower(orig.count - orig.ucount);
    if (bi == buckets.begin())
      lower.ucount = lower.count;
    else if (right == buckets.rbegin()->first)
      lower.ucount = lower.count/2;
    lower.stamp = orig.stamp;

    assert((use_sample || (bi != buckets.begin() &&
                           right != buckets.rbegin()->first)) &&
           "Splitting first or last bucket without a sample");

    T lright = use_sample ? sample :
                            right/2 + left/2; // distribute the division to avoid overflow.

    assert(value_classification<T>::finite(lright) &&
           "Not inserting finite lower bucket");

#ifdef DEBUG_HISTOGRAM
    std::cout << "inserting new bucket: " << lright << " right: " << right <<
            " left: " << left << " count: " << lower.count << "\n";
#endif

    std::pair<map_iterator, bool> r =
      buckets.insert(make_pair(lright, lower));
    if (!r.second) {
      // This is a duplicate bucket edge, conserve total counts by adding what
      // we would have inserted to the bucket already there (and mark all such
      // counts as upper counts).
      r.first->second.count += lower.count;
      r.first->second.ucount += lower.count;
      r.first->second.stamp = max(r.first->second.stamp, lower.stamp);
    }
  }

  map_iterator best_to_split(float &maxeps, map_iterator &bi, T sample) {
    maxeps = -std::numeric_limits<float>::infinity();
    map_iterator sbi = buckets.end();

    for (map_iterator i = buckets.begin(), b = buckets.begin(),
         e = buckets.end(); i != e; ++i) {
      if (!i->second.count)
        continue;

      T left = get_left_value(i);
      // Never split the first or last buckets (unless we have a sample).
      if (i != bi || sample == left) {
        if (i == b)
          continue;
        if (i->first == buckets.rbegin()->first)
          continue;
      }

      // Don't split a zero-width bin.
      if (i->first == left)
        continue;

      float avg = 0.5*i->second.count;
      float eps = fabs(i->second.ucount - avg) +
                  fabs((i->second.count - i->second.ucount) - avg);

      if (eps > 0 && eps > maxeps) {
        maxeps = eps;
        sbi = i;
      }
    }

    return sbi;
  }

  // Find the best bucket to merge with the bucket on the left.
  map_iterator best_to_merge(float &mineps) {
    mineps = std::numeric_limits<float>::infinity();
    map_iterator sbi = buckets.end();

    for (map_iterator i = buckets.begin(), b = buckets.begin(),
         e = buckets.end(); i != e; ++i) {
      // We're searching to the left, so skip the first bucket.
      if (i == b)
        continue;

      if (!i->second.count)
        continue;

      map_iterator l = i;
      --l;

      float avg = 0.25*i->second.count + 0.25*l->second.count;
      float eps = fabs(i->second.ucount - avg) + fabs((i->second.count - i->second.ucount) - avg) +
                  fabs(l->second.ucount - avg) + fabs((l->second.count - l->second.ucount) - avg);

      if (eps > 0 && eps < mineps) {
        mineps = eps;
        sbi = i;
      }
    }

    return sbi;
  }

protected:
  size_t nBuckets, nTotal, maxAge, curAge;
  T min_value, max_value;
  map_type buckets;
};

template <unsigned S, typename T>
struct uint_storage_ { typedef T type; };

template <typename T>
struct uint_storage_<1u, T>  { typedef unsigned char type; };
template <typename T>
struct uint_storage_<2u, T> { typedef uint16_t type; };
template <typename T>
struct uint_storage_<4u, T> { typedef uint32_t type; };
template <typename T>
struct uint_storage_<8u, T> { typedef uint64_t type; };

template <typename T>
struct uint_storage
  : uint_storage_<(std::numeric_limits<T>::is_signed ||
                   !std::numeric_limits<T>::is_integer)*sizeof(T), T> {};

template <unsigned S, typename T>
struct diff_ { typedef T type; };

template <typename T>
struct diff_<1u, T>  { typedef int16_t type; };
template <typename T>
struct diff_<2u, T> { typedef int32_t type; };
template <typename T>
struct diff_<4u, T> { typedef int64_t type; };

template <typename T>
struct diff
  : diff_<(!std::numeric_limits<T>::is_signed)*sizeof(T), T> {};

// The value bit predictor (using the previous value for context).
// This is the unique part of the algorithm.
template <typename T>
class numeric_value_predictor {
protected:
  typedef typename uint_storage<T>::type uint_type;
  typedef typename diff<T>::type diff_type;

public:
  numeric_value_predictor(size_t nBuckets)
    : h(nBuckets), has_prev(false), known(0), b(0) {}

  float predict() {
    float q = 0.5;

    if (has_prev && value_classification<T>::finite(prev)) {
      // Now we know the sign bit, but nothing else. This tells us whether (d < -data[i-1]),
      // and now we'd like to know what can be said about the first bit of the exponent.
      // The exponent stored, e, represents the exponent (e - 127), and so this first bit of
      // the exponent determines whether or not the overall exponent is positive or negative
      // (if this bit is set, then the exponent is at least one, if not, it is zero or
      // negative). If the bit is *not* set, then the number is:
      //   (-1)^sign (1 + [0, 1)) 2^0 or (-1)^sign 2 (meaning that the absolute value is at
      //   most 2).

      // largest and smallest numbers if the bit is not set.
      T large, small;
      // This is (31 - b) so that bit b remains clear.
      uint_type unknown_mask = (uint_type(1) << (sizeof(T)*8 - 1 - b)) - uint_type(1);
      uint_type ilarge = known | unknown_mask, ismall = known;
      memcpy(&large, &ilarge, sizeof(T));
      memcpy(&small, &ismall, sizeof(T));
	// Note: by testing for the unset bit, we can avoid creating a
	// fully-set exponent (which is a NaN).

      // IEEE floating point uses the all-set exponent to indicate NaN, Inf, etc.
      if (!value_classification<T>::finite(large)) {
        ilarge &= ~(uint_type(1) << (std::numeric_limits<T>::digits - 1));
        memcpy(&large, &ilarge, sizeof(T));
      }

      if (!value_classification<T>::finite(small)) {
        ismall &= ~(uint_type(1) << (std::numeric_limits<T>::digits - 1));
        memcpy(&small, &ismall, sizeof(T));
      }

      // data[i] = data[i-1] + d < max => d < max - data[i-1]
      // data[i] = data[i-1] + d > min => d > min - data[i-1]

      diff_type max_diff = diff_type(std::max(large, small)) - diff_type(prev);
      diff_type min_diff = diff_type(std::min(large, small)) - diff_type(prev);

      float q_unset = h.estimate_cdf(max_diff);
      if ((std::numeric_limits<diff_type>::is_integer && min_diff != std::numeric_limits<diff_type>::min()) ||
          (!std::numeric_limits<diff_type>::is_integer && min_diff != -std::numeric_limits<diff_type>::min())) {
        if (std::numeric_limits<diff_type>::is_integer)
          min_diff -= diff_type(1);
        else
          min_diff *= diff_type(1) - std::numeric_limits<diff_type>::epsilon();

        q_unset -= h.estimate_cdf(min_diff);
      }

#ifdef DEBUG_HISTOGRAM
      std::cout << "q_unset = " << q_unset << " (" << min_diff << ", " << max_diff << ")\n";
#endif

      uint_type set_mask = uint_type(1) << (sizeof(T)*8 - 1 - b);
      ilarge |= set_mask;
      ismall |= set_mask;
      memcpy(&large, &ilarge, sizeof(T));
      memcpy(&small, &ismall, sizeof(T));

      // IEEE floating point uses the all-set exponent to indicate NaN, Inf, etc.
      if (!value_classification<T>::finite(large)) {
        ilarge &= ~(uint_type(1) << (std::numeric_limits<T>::digits - 1));
        memcpy(&large, &ilarge, sizeof(T));
      }

      if (!value_classification<T>::finite(small)) {
        ismall &= ~(uint_type(1) << (std::numeric_limits<T>::digits - 1));
        memcpy(&small, &ismall, sizeof(T));
      }

      /*diff_type*/ max_diff = diff_type(std::max(large, small)) - diff_type(prev);
      /*diff_type*/ min_diff = diff_type(std::min(large, small)) - diff_type(prev);

      float q_set = h.estimate_cdf(max_diff);

      if ((std::numeric_limits<diff_type>::is_integer && min_diff != std::numeric_limits<diff_type>::min()) ||
          (!std::numeric_limits<diff_type>::is_integer && min_diff != -std::numeric_limits<diff_type>::min())) {
        if (std::numeric_limits<diff_type>::is_integer)
          min_diff -= diff_type(1);
        else
          min_diff *= diff_type(1) - std::numeric_limits<diff_type>::epsilon();

        q_set -= h.estimate_cdf(min_diff);
      }

      q = (q_set == 0 && q_unset == 0) ? 0.5f : q_set / (q_set + q_unset);

#ifdef DEBUG_HISTOGRAM
     std::cout << "q_set = " << q_set << " (" << min_diff << ", " << max_diff << ")\n";
     std::cout << "q = " << q << "\n";
#endif

      // Clamp the value; sometimes we end up slightly above or below zero.
      if (q > 1)
        q = 1;
      else if (q < 0)
        q = 0;
    }

    return q;
  }

  // Call set_bit after calling predict to record the actual bit value (either
  // the bit about to be encoded, or the bit just decoded).
  void set_bit(int y) {
    assert((y == 0 || y == 1) && "Invalid bit value");
    known |= uint_type(y) << (sizeof(T)*8 - 1 - b);

    // If we now have a completed known value, set it as the previous value,
    // and if we already had a previous value, then update the histogram.
    if (b == sizeof(T)*8 - 1) {
      T p;
      memcpy(&p, &known, sizeof(T));
      known = 0;

      if (has_prev && value_classification<T>::finite(p) &&
                      value_classification<T>::finite(prev)) {
        h.add(diff_type(p) - diff_type(prev));

#ifdef DEBUG_HISTOGRAM
        std::cout << "hist added: " << (diff_type(p) - diff_type(prev)) <<
                " (" << p << " vs. prev = " << prev << ")\n";
        h.print(std::cout);
        std::cout << "\n\n";
#endif
      }

      prev = p;
      has_prev = true;
    }

    if (++b == sizeof(T)*8)
      b = 0;
  }

  void print(std::ostream &os) const {
    h.print(os);
  }

protected:
#ifdef NO_SIMPLE_HISTOGRAM
  histogram<diff_type> h;
#else
  histogram_simple<diff_type> h;
#endif

  bool has_prev;
  T prev;

  uint_type known;
  int b;
};

template <typename T>
class numeric_value_predictor2 {
protected:
  typedef typename uint_storage<T>::type uint_type;

public:
  numeric_value_predictor2(size_t nBuckets)
    : h(nBuckets), known(0), b(0) {}

  float predict() {
    float q = 0.5;

    // Now we know the sign bit, but nothing else. This tells us whether (d < -data[i-1]),
    // and now we'd like to know what can be said about the first bit of the exponent.
    // The exponent stored, e, represents the exponent (e - 127), and so this first bit of
    // the exponent determines whether or not the overall exponent is positive or negative
    // (if this bit is set, then the exponent is at least one, if not, it is zero or
    // negative). If the bit is *not* set, then the number is:
    //   (-1)^sign (1 + [0, 1)) 2^0 or (-1)^sign 2 (meaning that the absolute value is at
    //   most 2).

    // largest and smallest numbers if the bit is not set.
    T large, small;
    // This is (31 - b) so that bit b remains clear.
    uint_type unknown_mask = (uint_type(1) << (sizeof(T)*8 - 1 - b)) - uint_type(1);
    uint_type ilarge = known | unknown_mask, ismall = known;
    memcpy(&large, &ilarge, sizeof(T));
    memcpy(&small, &ismall, sizeof(T));

    // IEEE floating point uses the all-set exponent to indicate NaN, Inf, etc.
    if (!value_classification<T>::finite(large)) {
      ilarge &= ~(uint_type(1) << (std::numeric_limits<T>::digits - 1));
      memcpy(&large, &ilarge, sizeof(T));
    }

    if (!value_classification<T>::finite(small)) {
      ismall &= ~(uint_type(1) << (std::numeric_limits<T>::digits - 1));
      memcpy(&small, &ismall, sizeof(T));
    }

    // data[i] = data[i-1] + d < max => d < max - data[i-1]
    // data[i] = data[i-1] + d > min => d > min - data[i-1]

    T max_diff = std::max(large, small);
    T min_diff = std::min(large, small);

    float q_unset = h.estimate_cdf(max_diff);

    if ((std::numeric_limits<T>::is_integer && min_diff != std::numeric_limits<T>::min()) ||
        (!std::numeric_limits<T>::is_integer && min_diff != -std::numeric_limits<T>::min())) {
      if (std::numeric_limits<T>::is_integer)
        min_diff -= T(1);
      else
        min_diff *= T(1) - std::numeric_limits<T>::epsilon();

      q_unset -= h.estimate_cdf(min_diff);
    }

#ifdef DEBUG_HISTOGRAM
    std::cout << "q_unset = " << q_unset << " (" << min_diff << ", " << max_diff << ")\n";
#endif

    uint_type set_mask = uint_type(1) << (sizeof(T)*8 - 1 - b);
    ilarge |= set_mask;
    ismall |= set_mask;
    memcpy(&large, &ilarge, sizeof(T));
    memcpy(&small, &ismall, sizeof(T));

    // IEEE floating point uses the all-set exponent to indicate NaN, Inf, etc.
    if (!value_classification<T>::finite(large)) {
      ilarge &= ~(uint_type(1) << (std::numeric_limits<T>::digits - 1));
      memcpy(&large, &ilarge, sizeof(T));
    }

    if (!value_classification<T>::finite(small)) {
      ismall &= ~(uint_type(1) << (std::numeric_limits<T>::digits - 1));
      memcpy(&small, &ismall, sizeof(T));
    }

    /*T*/ max_diff = std::max(large, small);
    /*T*/ min_diff = std::min(large, small);

    float q_set = h.estimate_cdf(max_diff);

    if ((std::numeric_limits<T>::is_integer && min_diff != std::numeric_limits<T>::min()) ||
        (!std::numeric_limits<T>::is_integer && min_diff != -std::numeric_limits<T>::min())) {
      if (std::numeric_limits<T>::is_integer)
        min_diff -= T(1);
      else
        min_diff *= T(1) - std::numeric_limits<T>::epsilon();

      q_set -= h.estimate_cdf(min_diff);
    }

    q = (q_set == 0 && q_unset == 0) ? 0.5f : q_set / (q_set + q_unset);

#ifdef DEBUG_HISTOGRAM
   std::cout << "q_set = " << q_set << " (" << min_diff << ", " << max_diff << ")\n";
   std::cout << "q = " << q << "\n";
#endif

    // Clamp the value; sometimes we end up slightly above or below zero.
    if (q > 1)
      q = 1;
    else if (q < 0)
      q = 0;

    return q;
  }

  // Call set_bit after calling predict to record the actual bit value (either
  // the bit about to be encoded, or the bit just decoded).
  void set_bit(int y) {
    assert((y == 0 || y == 1) && "Invalid bit value");
    known |= uint_type(y) << (sizeof(T)*8 - 1 - b);

    // If we now have a completed known value, set it as the previous value,
    // and if we already had a previous value, then update the histogram.
    if (b == sizeof(T)*8 - 1) {
      T p;
      memcpy(&p, &known, sizeof(T));
      known = 0;

      if (value_classification<T>::finite(p)) {
        h.add(p);

#ifdef DEBUG_HISTOGRAM
        std::cout << "hist added: " << p << "\n";
        h.print(std::cout);
        std::cout << "\n\n";
#endif
      }
    }

    if (++b == sizeof(T)*8)
      b = 0;
  }

  void print(std::ostream &os) const {
    h.print(os);
  }

protected:
#ifdef NO_SIMPLE_HISTOGRAM
  histogram<T> h;
#else
  histogram_simple<T> h;
#endif

  uint_type known;
  int b;
};

// Logistic model mixing; see:
// http://en.wikipedia.org/wiki/Context_mixing
template <typename T>
class combination_predictor {
protected:
  static const int nModels = 2;

public:
  combination_predictor(size_t nBuckets, float e = 0.0015/3)
    : pred(nBuckets), pred2(nBuckets), eta(e) {
    for (int i = 0; i < nModels; ++i)
      for (b = 0; b < (int) sizeof(T)*8; ++b)
        w[i][b] = 0;

    b = 0;
  }

  float predict() {
    float q1 = stretch(pred.predict()),
          q2 = stretch(pred2.predict());

    p[0][b] = q1;
    p[1][b] = q2;

    p1 = squash(w[0][b]*q1 + w[1][b]*q2);
    return p1;
  }

  void set_bit(int y) {
    pred.set_bit(y);
    pred2.set_bit(y);

    for (int i = 0; i < nModels; ++i)
      w[i][b] += eta * p[i][b] * (y - p1);

    if (++b == sizeof(T)*8)
      b = 0;
  }

  void print(std::ostream &os) const {
    os << "Model 1:\n";
    for (int c = 0; c < (int) sizeof(T)*8; ++c)
      os << "  bit " << c << ": " << w[0][c] << "\n";
    os << "\n";
    pred.print(os);

    os << "Model 2:\n";
    for (int c = 0; c < (int) sizeof(T)*8; ++c)
      os << "  bit " << c << ": " << w[1][c] << "\n";
    os << "\n";
    pred2.print(os);
  }

protected:
  float stretch(float x) {
    if (x == 1.0f)
      x -= std::numeric_limits<float>::epsilon();
    else if (x == 0)
      x += std::numeric_limits<float>::min();

    return log(double(x / (1.0f - x)));
  }

  float squash(float x) {
    return 1.0f / (1.0f + exp(double(-x)));
  }

protected:
  numeric_value_predictor<T>  pred;
  numeric_value_predictor2<T> pred2;

  int b;
  float eta;
  float w[nModels][sizeof(T)*8];
  float p[nModels][sizeof(T)*8];
  float p1;
};

// The compresor for slowly-varying data.
template <typename T>
struct sld_compressor {
protected:
  typedef typename uint_storage<T>::type uint_type;

public:
  // Extract bit b (where b = 0 is the high-order bit).
  static int extract_bit(T f, int b) {
    uint_type v;
    memcpy(&v, &f, sizeof(T));
    return (v >> (sizeof(T)*8 - 1 - b)) & 0x1;
  }

  // Extract bit b (where b = 0 is the high-order bit).
  static void set_bit(T &f, int b, int y) {
    assert((y == 0 || y == 1) && "Invalid bit value");
    uint_type v;
    memcpy(&v, &f, sizeof(T));
    v |= (uint_type(y) << (sizeof(T)*8 - 1 - b));
    memcpy(&f, &v, sizeof(T));
  }

	template <typename IT>
	void decompressData(IT data_begin, IT data_end, size_t nBuckets, std::vector<unsigned char> &cdata);

	template <typename IT>
	void compressData(IT data_begin, IT data_end, size_t nBuckets, std::vector<unsigned char> &cdata);

	template <typename IT>
	void compressData_omp(IT data_begin, IT data_end, size_t nBuckets, size_t blockSize, std::vector<unsigned char> &cdata);

	template <typename IT>
	void decompressData_omp(IT data_begin, IT data_end, size_t nBuckets, size_t blockSize, std::vector<unsigned char> &cdata, size_t cdstart);
};


template <typename T> template <typename IT>
void sld_compressor<T>::decompressData(IT data_begin, IT data_end, size_t nBuckets, std::vector<unsigned char> &cdata)
{
	combination_predictor<T> pred(nBuckets);
	arithmetic_coding ac(cdata, false);

	for (IT data = data_begin; data != data_end; ++data)
	{
		T v = 0;
		for (int b = 0; b < (int) sizeof(T)*8; ++b)
		{
			float q = pred.predict();
			int y = ac.decode(q);

			set_bit(v, b, y);
			pred.set_bit(y);
		}

		*data = v;
	}
}


template <typename T> template <typename IT>
void sld_compressor<T>::compressData(IT data_begin, IT data_end, size_t nBuckets, std::vector<unsigned char> &cdata)
{
	combination_predictor<T> pred(nBuckets);
	arithmetic_coding ac(cdata);

	for (IT data = data_begin; data != data_end; ++data)
	{
		for (int b = 0; b < (int) sizeof(T)*8; ++b)
		{
			float q = pred.predict();
			int y = extract_bit(*data, b);

			ac.encode(q, y);
			pred.set_bit(y);
		}
	}
	ac.flush();
}


template <typename T> template <typename IT>
void sld_compressor<T>::compressData_omp(IT data_begin, IT data_end, size_t nBuckets, size_t blockSize, std::vector<unsigned char> &cdata)
{
	size_t Np = std::distance(data_begin, data_end);
	size_t nBlocks = (Np + blockSize - 1)/blockSize;

	size_t BSAoff = cdata.size();
	cdata.resize(cdata.size() + sizeof(uint64_t)*nBlocks);

	std::vector<std::vector<unsigned char> > cdb(nBlocks);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
	for (int i = 0; i < (int) nBlocks; ++i)
	{
		IT ldata_begin = data_begin, ldata_end = data_begin;
		std::advance(ldata_begin, i*blockSize);
		if (i == (int) nBlocks - 1)
			ldata_end = data_end;
		else
			std::advance(ldata_end, (i+1)*blockSize);

		compressData(ldata_begin, ldata_end, nBuckets, cdb[i]);
	}

	for (int i = 0; i < (int) nBlocks; ++i)
	{
		size_t off = cdata.size();
		((uint64_t*) &cdata[BSAoff])[i] = off;
		cdata.resize(cdata.size() + cdb[i].size());
		std::copy(cdb[i].begin(), cdb[i].end(), cdata.begin() + off);
	}
}


template <typename T> template <typename IT>
void sld_compressor<T>::decompressData_omp(IT data_begin, IT data_end, size_t nBuckets, size_t blockSize, std::vector<unsigned char> &cdata, size_t cdstart)
{
	size_t Np = std::distance(data_begin, data_end);
	size_t nBlocks = (Np + blockSize - 1)/blockSize;

	uint64_t *BSA = (uint64_t *) &cdata[cdstart];

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
	for (int i = 0; i < (int) nBlocks; ++i)
	{
		IT ldata_begin = data_begin, ldata_end = data_begin;
		std::advance(ldata_begin, i*blockSize);
		if (i == (int) nBlocks - 1)
			ldata_end = data_end;
		else
			std::advance(ldata_end, (i+1)*blockSize);

		std::vector<unsigned char> cdb(cdata.begin() + BSA[i], i+1 < (int) nBlocks ? cdata.begin() + BSA[i+1] : cdata.end());
		decompressData(ldata_begin, ldata_end, nBuckets, cdb);
	}
}

}

#endif // GENERICIOCOMPRESS_H

