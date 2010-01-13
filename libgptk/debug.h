#ifdef DEBUG // DEBUG MODE ENABLED                 
 	#define debug_msg(msg) cout << "** DEBUG: " << msg << endl;
#else          // DEBUG MODE DISABLED                        
 	#define debug_msg(msg)
 	#ifndef NDEBUG
  		#define NDEBUG              // Disable assertions
 	#endif
#endif

#include <cassert>
