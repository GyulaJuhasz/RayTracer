#ifndef PROGRAM_H
#define PROGRAM_H

void initialize();

void onDisplay();

void onResize(int w, int h);

void onMouseEvent(int button, int state, int x, int y);

void onKeyboardEvent(unsigned char key, int x, int y);

#endif // !PROGRAM_H