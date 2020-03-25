#include "GaudiKernel/DeclareFactoryEntries.h"
#include "TestpiAlg/Testpi.h"

DECLARE_ALGORITHM_FACTORY( Testpi )

DECLARE_FACTORY_ENTRIES( TestpiAlg ) {
  DECLARE_ALGORITHM(Testpi);
}

