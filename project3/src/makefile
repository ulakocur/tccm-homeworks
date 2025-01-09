CC = gcc

SOURCES = main.c functions.c
TARGET = md_simulation

$(TARGET): $(SOURCES)
	$(CC)  -o $@ $^ 
	mv $(TARGET) ../
	
clean:
	rm -f $(TARGET)

