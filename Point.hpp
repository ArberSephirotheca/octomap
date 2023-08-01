    // 3D Point structure
    struct Point
    {
        double x, y, z;
        Point(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
        Point(): x(0.0), y(0.0), z(0.0) {}
    };