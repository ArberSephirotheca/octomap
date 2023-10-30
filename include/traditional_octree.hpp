#include "traditional_octree_base.hpp"
using namespace std;


// Octree class
class OctreeNode : public OctreeNodeBase
{

    public:
    Point* point;
 
    // Represent the boundary of the cube
    OctreeNode* children[8];
    Point *topLeftFront, *bottomRightBack;
    OctreeNode()
    {
        point = new Point();
    }
    // Constructor with three arguments
    OctreeNode(float x, float y, float z)
    {
        // To declare point node
        point = new Point(x, y, z);
    }
    OctreeNode(float x1, float y1, float z1, float x2, float y2, float z2)
    {
        // This use to construct Octree
        // with boundaries defined
        if (x2 < x1 || y2 < y1 || z2 < z1)
        {
            cout << "boundary points are not valid" << endl;
            return;
        }
        point = nullptr;
        topLeftFront = new Point(x1, y1, z1);
        bottomRightBack = new Point(x2, y2, z2);
        for(int i = 0; i < 8; ++i){
            children[i] = nullptr;
        }
        for (int i = TopLeftFront;
             i <= BottomLeftBack;
             ++i)
            children[i] = new OctreeNode();
    }
    // Function that returns true if the point
    // (x, y, z) exists in the octree
     bool find(int x, int y, int z)
    {
        // If point is out of bound
        if (x < topLeftFront->x
            || x > bottomRightBack->x
            || y < topLeftFront->y
            || y > bottomRightBack->y
            || z < topLeftFront->z
            || z > bottomRightBack->z)
            return 0;
 
        // Otherwise perform binary search
        // for each ordinate
        int midx = (topLeftFront->x
                    + bottomRightBack->x)
                   / 2;
        int midy = (topLeftFront->y
                    + bottomRightBack->y)
                   / 2;
        int midz = (topLeftFront->z
                    + bottomRightBack->z)
                   / 2;
 
        int pos = -1;
 
        // Deciding the position
        // where to move
        if (x <= midx) {
            if (y <= midy) {
                if (z <= midz)
                    pos = TopLeftFront;
                else
                    pos = TopLeftBottom;
            }
            else {
                if (z <= midz)
                    pos = BottomLeftFront;
                else
                    pos = BottomLeftBack;
            }
        }
        else {
            if (y <= midy) {
                if (z <= midz)
                    pos = TopRightFront;
                else
                    pos = TopRightBottom;
            }
            else {
                if (z <= midz)
                    pos = BottomRightFront;
                else
                    pos = BottomRightBack;
            }
        }
 
        // If an internal node is encountered
        if (children[pos]->point == nullptr) {
            return children[pos]->find(x, y, z);
        }
        // If an empty node is encountered
        if (children[pos]->point->x == -1.0) {
            return false;
        }
        else { 
            // If node is found with
            // the given value
            if (x == children[pos]->point->x
                && y == children[pos]->point->y
                && z == children[pos]->point->z)
                return 1;
        }
        return 0;
    }

    void clear()
    {
        for (int i = 0; i < 8; ++i)
        {
            if (children[i])
            {
                children[i]->clear();
                delete children[i];
            }
        }
        delete point;
        delete topLeftFront;
        delete bottomRightBack;
    }
};

class Octree : public OctreeBase
{
    OctreeNode *_root;

public:
    Octree(float x1, float y1, float z1, float x2, float y2, float z2)
    {
        _root = new OctreeNode(x1, y1, z1, x2, y2, z2);
    }
    ~Octree()
    {
        if (_root)
        {
            _root->clear();
        }
        delete _root;
    }
    void insert(float x, float y, float z)
    {
        // If the point already exists in the octree
        if (find(x, y, z) == true)
        {
            return;
        }
        insert(_root, x, y, z);
    }
    bool find(float x, float y, float z){
        return _root->find(x, y, z);
    }
    void insert(OctreeNode *node, float x, float y, float z){
        // If the point is out of bounds
        if (x < node->topLeftFront->x || x > node->bottomRightBack->x || y < node->topLeftFront->y || y > node->bottomRightBack->y || z < node->topLeftFront->z || z > node->bottomRightBack->z)
        {
            return;
        }

        // Binary search to insert the point
        float midx = (node->topLeftFront->x + node->bottomRightBack->x) / 2;
        float midy = (node->topLeftFront->y + node->bottomRightBack->y) / 2;
        float midz = (node->topLeftFront->z + node->bottomRightBack->z) / 2;

        int pos = -1;

        // Checking the octant of
        // the point
        if (x <= midx)
        {
            if (y <= midy)
            {
                if (z <= midz)
                    pos = TopLeftFront;
                else
                    pos = TopLeftBottom;
            }
            else
            {
                if (z <= midz)
                    pos = BottomLeftFront;
                else
                    pos = BottomLeftBack;
            }
        }
        else
        {
            if (y <= midy)
            {
                if (z <= midz)
                    pos = TopRightFront;
                else
                    pos = TopRightBottom;
            }
            else
            {
                if (z <= midz)
                    pos = BottomRightFront;
                else
                    pos = BottomRightBack;
            }
        }

        // If an internal node is encountered
        if (node->children[pos]->point == nullptr)
        {
            //std::cout<<"internal point"<<std::endl;
            insert(node->children[pos], x, y, z);
            return;
        }

        // If an empty node is encountered
        else if (node->children[pos]->point->x == -1)
        {
            //std::cout <<"empty point"<<std::endl;
            delete node->children[pos];
            node->children[pos] = new OctreeNode(x, y, z);
            return;
        }
        else
        {
            float x_ = node->children[pos]->point->x,
                  y_ = node->children[pos]->point->y,
                  z_ = node->children[pos]->point->z;
            delete node->children[pos];
            node->children[pos] = nullptr;
            if (pos == TopLeftFront)
            {
                node->children[pos] = new OctreeNode(node->topLeftFront->x,
                                               node->topLeftFront->y,
                                               node->topLeftFront->z,
                                               midx,
                                               midy,
                                               midz);
            }

            else if (pos == TopRightFront)
            {
                node->children[pos] = new OctreeNode(midx + 1,
                                               node->topLeftFront->y,
                                               node->topLeftFront->z,
                                               node->bottomRightBack->x,
                                               midy,
                                               midz);
            }
            else if (pos == BottomRightFront)
            {
                node->children[pos] = new OctreeNode(midx + 1,
                                               midy + 1,
                                               node->topLeftFront->z,
                                               node->bottomRightBack->x,
                                               node->bottomRightBack->y,
                                               midz);
            }
            else if (pos == BottomLeftFront)
            {
                node->children[pos] = new OctreeNode(node->topLeftFront->x,
                                               midy + 1,
                                               node->topLeftFront->z,
                                               midx,
                                               node->bottomRightBack->y,
                                               midz);
            }
            else if (pos == TopLeftBottom)
            {
                node->children[pos] = new OctreeNode(node->topLeftFront->x,
                                               node->topLeftFront->y,
                                               midz + 1,
                                               midx,
                                               midy,
                                               node->bottomRightBack->z);
            }
            else if (pos == TopRightBottom)
            {
                node->children[pos] = new OctreeNode(midx + 1,
                                               node->topLeftFront->y,
                                               midz + 1,
                                               node->bottomRightBack->x,
                                               midy,
                                               node->bottomRightBack->z);
            }
            else if (pos == BottomRightBack)
            {
                node->children[pos] = new OctreeNode(midx + 1,
                                               midy + 1,
                                               midz + 1,
                                               node->bottomRightBack->x,
                                               node->bottomRightBack->y,
                                               node->bottomRightBack->z);
            }
            else if (pos == BottomLeftBack)
            {
                node->children[pos] = new OctreeNode(node->topLeftFront->x,
                                               midy + 1,
                                               midz + 1,
                                               midx,
                                               node->bottomRightBack->y,
                                               node->bottomRightBack->z);
            }
            insert(node->children[pos], x_, y_, z_);
            insert(node->children[pos], x, y, z);
        }
    }
    void clear()
    {
        for (int i = 0; i < 8; ++i)
        {
            if (_root->children[i])
            {
                _root->children[i]->clear();
                delete _root->children[i];
            }
        }
    }
};