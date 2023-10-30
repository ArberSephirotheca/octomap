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
    float x;
    float y;
    float z;
    Point()
        : x(-1.0), y(-1.0), z(-1.0)
    {
    }
 
    Point(float a, float b, float c)
        : x(a), y(b), z(c)
    {
    }
};

class OctreeNodeBase{
};

class OctreeBase{
    public:
    virtual ~OctreeBase(){};
    virtual void insert(float x,
                float y,
                float z){};
    virtual bool find(float x, float y, float z){};
    virtual void clear() {};
};