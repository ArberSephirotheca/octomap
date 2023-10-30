// Structure of a point
#pragma once

#define TopLeftFront 0
#define TopRightFront 1
#define BottomRightFront 2
#define BottomLeftFront 3
#define TopLeftBottom 4
#define TopRightBottom 5
#define BottomRightBack 6
#define BottomLeftBack 7

struct Point {
    int x;
    int y;
    int z;
    Point()
        : x(-1), y(-1), z(-1)
    {
    }
 
    Point(int a, int b, int c)
        : x(a), y(b), z(c)
    {
    }
};

class OctreeNodeBase{
};

class OctreeBase{
    public:
    virtual ~OctreeBase(){};
    virtual void insert(int x,
                int y,
                int z){};
    virtual bool find(int x, int y, int z){};
    virtual void clear() {};
};