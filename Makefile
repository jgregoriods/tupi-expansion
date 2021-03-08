TARGET = expand

CC = g++ -std=c++17
CFLAGS = -Ofast
LINKER = g++ -std=c++17

SRCDIR = src
OBJDIR = obj

SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

$(TARGET): $(OBJECTS)
	$(LINKER) -o $(TARGET) $(OBJECTS)

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	@$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm $(OBJECTS)
	rm expand
