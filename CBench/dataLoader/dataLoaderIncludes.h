#include "HACC/HACCDataLoader.hpp"

#ifdef CBENCH_HAS_BINARY
	#include "Binary/BinaryDataLoader.hpp"
#endif

#ifdef CBENCH_HAS_NYX
	#include "NYX/NYXDataLoader.hpp"
#endif

#ifdef CBENCH_HAS_VTK
	#include "VTK/VTKDataLoader.hpp"
#endif

#ifdef CBENCH_HAS_GDA
	#include "VPIC_GDA/GDADataLoader.hpp"
#endif

#ifdef CBENCH_HAS_AMReX
	#include "AMReX/PltDataLoader.hpp"
#endif
