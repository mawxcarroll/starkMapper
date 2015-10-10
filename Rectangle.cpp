#include <iostream>
#include "Rectangle.h"
using namespace std;

Rectangle::Rectangle(int w, int h)
{
    width = w;
    height = h;
}

int Rectangle::area()
{
    return width*height;
}

void Rectangle::set_values (int x, int y)
{
    width = x;
    height = y;
}
