class Code
{
public:
    Code() : code_(0), depth_(0) {}
    Code(uint32_t code, int depth = 0) : code_(code), depth_(depth) {
    }

    uint32_t getCode() const { return code_; }
    // Get all children that this code can have from this code's depth to depth 0
    std::vector<Code> getAllChildren() const
    {
        std::vector<Code> children;
        uint32_t max = 8 << (3 * depth_);
        for (uint32_t i = 0; i < max; ++i)
        {
            children.emplace_back(code_ + i, 0);
        }
        return children;
    }

    int getChildIndex(int depth) const
    {
        return (code_ >> (3 * (depth - 1))) & 0x7;
    }

    Code getChild(int index)
    {
        if (0 == depth_)
        {
            return *this;
        }

        int child_depth = depth_ - 1;
        return Code(code_ + (static_cast<uint32_t>(index) << static_cast<uint32_t>(3 * child_depth)), child_depth);
    }

    std::vector<Code> getChildren() const
    {
        std::vector<Code> children;
        if (0 == depth_)
        {
            return children;
        }

        int child_depth = depth_ - 1;
        uint32_t offset = 3 * child_depth;
        for (uint32_t i = 0; i < 8; ++i)
        {
            children.emplace_back(code_ + (i << offset), child_depth);
        }
        return children;
    }

    int getDepth() const
    {
        return depth_;
    }

    // see https://developer.nvidia.com/blog/thinking-parallel-part-iii-tree-construction-gpu/
    unsigned int expandBits(unsigned int v)
    {
        v = (v * 0x00010001u) & 0xFF0000FFu;
        v = (v * 0x00000101u) & 0x0F00F00Fu;
        v = (v * 0x00000011u) & 0xC30C30C3u;
        v = (v * 0x00000005u) & 0x49249249u;
        return v;
    }

    // Calculates a 30-bit Morton code for the
    // given 3D point located within the unit cube [0,1].
    unsigned int morton3D(const Point &point)
    {
        double x = point.x;
        double y = point.y;
        double z = point.z;
        x = std::min(std::max(x * 1024.0, 0.0), 1023.0);
        y = std::min(std::max(y * 1024.0, 0.0), 1023.0);
        z = std::min(std::max(z * 1024.0, 0.0), 1023.0);
        unsigned int xx = expandBits((unsigned int)x);
        unsigned int yy = expandBits((unsigned int)y);
        unsigned int zz = expandBits((unsigned int)z);
        return xx * 4 + yy * 2 + zz;
    }
    // problem: maybt change it to uint64_t?
    private:
    uint32_t code_;
    int depth_;
};

void radixSortMorton(std::vector<uint32_t>& mortonCodes) {
    std::vector<unsigned int> temp(mortonCodes.size());

    // Perform radix sort on each bit position
    for (int shift = 0; shift < 30; ++shift) {
        std::vector<int> buckets[2];
        for (unsigned int code : mortonCodes) {
            int bucketIndex = (code >> shift) & 1;
            buckets[bucketIndex].push_back(code);
        }
        mortonCodes.clear();
        for (int i = 0; i < 2; ++i) {
            mortonCodes.insert(mortonCodes.end(), buckets[i].begin(), buckets[i].end());
        }
    }
}