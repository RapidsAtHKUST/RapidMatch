set(ND_SRC
        nd.h
        bucket.h
        bucket.cpp
        hierarchy.cpp
        k3-4.cpp
        kcore.cpp
        ktruss.cpp
        util.cpp
        nd_interface.cpp
        nd_interface.h)

add_library(nd SHARED ${ND_SRC})