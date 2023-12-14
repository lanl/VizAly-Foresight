#include "HACC/HACCDataLoader.hpp"

#ifdef CBENCH_HAS_NYX
	#include "NYX/NYXDataLoader.hpp"
#endif

#ifdef CBENCH_HAS_VTK
	#include "VTK/VTKDataLoader.hpp"
	#include "VTP/VTPDataLoader.hpp"
#endif

#ifdef CBENCH_HAS_GDA
	#include "VPIC_GDA/GDADataLoader.hpp"
#endif

#ifdef CBENCH_HAS_RAW
	#include "Generic_Binary/GenericBinaryLoader.hpp"
#endif