#include <type_traits>
namespace redwood
{
    template <typename T>
    void trash(T &&t)
    {
        static_assert(!std::is_same<T, void>::value, "");
    }
    class PID
    {
    public:
        static int get_pid();
        static int get_parent_pid();
    };
}