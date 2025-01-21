TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        elem/Element.cpp \
        elem/ElementQuad2D.cpp \
        elem/ElementQuad3D.cpp \
        elem/ShapeQuad2D.cpp \
        elem/ShapeQuad3D.cpp \
        elem/StressContainer.cpp \
        fea/FE.cpp \
        main.cpp \
        material/ElasticMaterial.cpp \
        material/Material.cpp \
        model/Dof.cpp \
        model/ElemFaceLoad.cpp \
        model/FeLoadData.cpp \
        model/FeModel.cpp \
        model/FeModelData.cpp \
        model/FeStress.cpp \
        solver/Solver.cpp \
        solver/SolverLDU.cpp \
        solver/SolverPCG.cpp \
        util/FePrintWriter.cpp \
        util/GaussRule.cpp \
        util/UTIL.cpp

HEADERS += \
    elem/Element.h \
    elem/ElementQuad2D.h \
    elem/ElementQuad3D.h \
    elem/ShapeQuad2D.h \
    elem/ShapeQuad3D.h \
    elem/StressContainer.h \
    fea/FE.h \
    material/ElasticMaterial.h \
    material/Material.h \
    model/Dof.h \
    model/ElemFaceLoad.h \
    model/FeLoadData.h \
    model/FeModel.h \
    model/FeModelData.h \
    model/FeStress.h \
    solver/Solver.h \
    solver/SolverLDU.h \
    solver/SolverPCG.h \
    util/FePrintWriter.h \
    util/GaussRule.h \
    util/UTIL.h
