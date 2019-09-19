#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-MacOSX
CND_DLIB_EXT=dylib
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/1946826759/Comp.o \
	${OBJECTDIR}/_ext/1946826759/Romb.o \
	${OBJECTDIR}/_ext/1946826759/algorithm.o \
	${OBJECTDIR}/_ext/1946826759/main.o \
	${OBJECTDIR}/_ext/1946826759/newComp.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/simpson_revised

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/simpson_revised: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/simpson_revised ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/_ext/1946826759/Comp.o: /Users/masha/NetBeansProjects/Simpson_revised/src/Comp.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1946826759
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1946826759/Comp.o /Users/masha/NetBeansProjects/Simpson_revised/src/Comp.cc

${OBJECTDIR}/_ext/1946826759/Romb.o: /Users/masha/NetBeansProjects/Simpson_revised/src/Romb.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1946826759
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1946826759/Romb.o /Users/masha/NetBeansProjects/Simpson_revised/src/Romb.cpp

${OBJECTDIR}/_ext/1946826759/algorithm.o: /Users/masha/NetBeansProjects/Simpson_revised/src/algorithm.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1946826759
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1946826759/algorithm.o /Users/masha/NetBeansProjects/Simpson_revised/src/algorithm.cpp

${OBJECTDIR}/_ext/1946826759/main.o: /Users/masha/NetBeansProjects/Simpson_revised/src/main.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1946826759
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1946826759/main.o /Users/masha/NetBeansProjects/Simpson_revised/src/main.cc

${OBJECTDIR}/_ext/1946826759/newComp.o: /Users/masha/NetBeansProjects/Simpson_revised/src/newComp.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1946826759
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1946826759/newComp.o /Users/masha/NetBeansProjects/Simpson_revised/src/newComp.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/simpson_revised

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
