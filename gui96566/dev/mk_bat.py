#!/usr/bin/python

import os
import shutil

path = "./crypto_sign/"
tar1 = [ "quartz10312934/" , "gui95955/" , "gui941744/" , "gui127946/" ]
cfg1 = [ "#define QUARTZ103\n" , "#define QUARTZ95\n" , "#define QUARTZ94\n" , "#define QUARTZ127\n"  ]
pmod1 = [ "ref" , "pclmulqdq" , "haswell" ]

tar2 = [ "gui96566/" , "gui128946/" ]
cfg2 = [ "#define QUARTZ96\n" , "#define QUARTZ128\n" ]
pmod2 = [ "ref" , "pclmulqdq" , "haswell" , "pshufb" , "neon" ]

whitelist0 = [ "designers" , "implementors" , "description" ]
whitelist1 = [ "api.h" , "blas.h" , "blas.cpp" , "config.h" , "gf2ext.hpp" , "gf2ext.cpp" , "gf2ext_poly.hpp" ,
"gf2ext_poly.cpp" , "mpkc.hpp" , "poly_lib.h" , "quartz_core.h" , "quartz_core.hpp" , "quartz.hpp" , "quartz.h" ,
"quartz.cpp"  , "sizes.h" , "usehash.cpp" , "crand.h" ]
#, "crypto_sign.h" , "crypto_hash_sha256.h" , "tester-test.cpp" ]

whitelist2 = [ "gf2ext-iso.h" , "gf2ext-iso.cpp" , "gf2ext-sse.hpp" , "gf2ext-sse.cpp" , "blas-sse.h" ]
whitelist3 = [ "gf2ext-iso.h" , "gf2ext-iso.cpp" , "gf2ext-neon.hpp" , "gf2ext-neon.cpp" ]


for p in tar1:
  for q in pmod1:
	dest = path+p+q
	if not os.path.exists(dest): os.makedirs(dest)

for p in tar2:
  for q in pmod2:
	dest = path+p+q
	if not os.path.exists(dest): os.makedirs(dest)


#src_files = os.listdir("./")
#for file_name in src_files:
#    jump = False;
#    for i in ignore:
#      if( -1 != file_name.find(i) ) : jump = True;
#    if( jump ) : continue;

for p in range(0,len(tar1)):
    for q in pmod1:
        tt = tar1[p]
	dest = path+tt+q
        f = open(dest + "/scheme.h" , 'w')
        f.write( cfg1[p] )
        f.close()
        f = open(dest + "/run_config.h" , 'w')
        f.write( "#ifndef _RUN_CONFIG_H_\n#define _RUN_CONFIG_H_\n\n" )
        if ( "pclmulqdq" == q ) :
           f.write("#define CONFIG_SSE_VEC\n")
           f.write("#define CONFIG_HAS_PCLMULQDQ\n")
           for file_name in whitelist2:
              full_file_name = os.path.join("./", file_name)
              shutil.copy(full_file_name, path+tt+q+"/")
        if ( "haswell" == q ) :
           f.write("#define CONFIG_SSE_VEC\n")
           f.write("#define CONFIG_HAS_PCLMULQDQ\n")
           f.write("#define CONFIG_FAST_PCLMULQDQ\n")
           for file_name in whitelist2:
              full_file_name = os.path.join("./", file_name)
              shutil.copy(full_file_name, path+tt+q+"/")
        if ( "pshufb" == q ) :
           f.write("#define CONFIG_SSE_VEC\n")
           f.write("#define CONFIG_PSHUFB\n")
           for file_name in whitelist2:
              full_file_name = os.path.join("./", file_name)
              shutil.copy(full_file_name, path+tt+q+"/")
        if ( "neon" == q ) :
           f.write("#define CONFIG_NEON\n")
           for file_name in whitelist3:
              full_file_name = os.path.join("./", file_name)
              shutil.copy(full_file_name, path+tt+q+"/")
        f.write( "\n\n#endif\n" )
        f.close()
        for file_name in whitelist1:
            full_file_name = os.path.join("./", file_name)
            shutil.copy(full_file_name, path+tt+q+"/")
        for file_name in whitelist0:
            full_file_name = os.path.join("./", file_name)
            shutil.copy(full_file_name, path+tt+q+"/")



for p in range(0,len(tar2)):
    for q in pmod2:
        tt = tar2[p]
	dest = path+tt+q
        f = open(dest + "/scheme.h" , 'w')
        f.write( cfg2[p] )
        f.close()
        f = open(dest + "/run_config.h" , 'w')
        f.write( "#ifndef _RUN_CONFIG_H_\n#define _RUN_CONFIG_H_\n\n" )
        if ( "pclmulqdq" == q ) :
           f.write("#define CONFIG_SSE_VEC\n")
           f.write("#define CONFIG_HAS_PCLMULQDQ\n")
           for file_name in whitelist2:
              full_file_name = os.path.join("./", file_name)
              shutil.copy(full_file_name, path+tt+q+"/")
        if ( "haswell" == q ) :
           f.write("#define CONFIG_SSE_VEC\n")
           f.write("#define CONFIG_HAS_PCLMULQDQ\n")
           f.write("#define CONFIG_FAST_PCLMULQDQ\n")
           for file_name in whitelist2:
              full_file_name = os.path.join("./", file_name)
              shutil.copy(full_file_name, path+tt+q+"/")
        if ( "pshufb" == q ) :
           f.write("#define CONFIG_SSE_VEC\n")
           f.write("#define CONFIG_PSHUFB\n")
           for file_name in whitelist2:
              full_file_name = os.path.join("./", file_name)
              shutil.copy(full_file_name, path+tt+q+"/")
        if ( "neon" == q ) :
           f.write("#define CONFIG_NEON\n")
           for file_name in whitelist3:
              full_file_name = os.path.join("./", file_name)
              shutil.copy(full_file_name, path+tt+q+"/")
        f.write( "\n\n#endif\n" )
        f.close()
        for file_name in whitelist1:
            full_file_name = os.path.join("./", file_name)
            shutil.copy(full_file_name, path+tt+q+"/")
        for file_name in whitelist0:
            full_file_name = os.path.join("./", file_name)
            shutil.copy(full_file_name, path+tt+q+"/")

