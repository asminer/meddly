
# $Id$

#
# Right now - proof of concept makefile
#

all: thing.jar

SRC= \
      info/LeafInfo.java      \
      info/ForestInfo.java    \
      logging/Logger.java     \
      logic/ForestInfoParser.java     \
      main/ApplicationExecution.java



OBJ= $(SRC:%.java=bin/%.class)

thing.jar: $(OBJ)
	jar cf thing.jar $(OBJ)

#
# Generic rule
#

bin/%.class: src/%.java
	javac -classpath bin/ -d bin/ $<


#
# Dependencies
#

bin/logic/ForestInforParser.class:  bin/info/LeafInfo.class \
                                    bin/info/ForestInfo.class

bin/main/ApplicationExecution.class:  bin/info/LeafInfo.class \
                                      bin/logic/ForestInforParser.class
