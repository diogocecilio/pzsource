// Emacs will be in -*- Mode: c++ -*-
//
// ***************** DO NOT REMOVE THIS BANNER *****************
//
// SUMMARY: Tools for Automatic Differentiaton (order 1)
// RELEASE: 0.1     
// USAGE  : You may copy freely these files and use it for    
//          teaching or research. These or part of these may   
//          not be sold or used for a commercial purpose with- 
//          out our consent : fax (33)1 44 27 72 00
//			
// AUTHOR : Nicolas Di cesare
// ORG    :          
// E-MAIL : Nicolas.Dicesare@ann.jussieu.fr     
//
// ORIG-DATE: September 97
// LAST-MOD : 20/12/97
// ************************************************************     
#ifndef _error_h
#define _error_h    

#include <cstdlib>
#include <iostream>
#include <string>


inline void error(const std::string &msg)
{
  if(msg.size()) std::cout << msg.c_str() << std::endl;
}

     
#endif
