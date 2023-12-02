
#include <iostream>

class myindex {
    public:
        myindex(int i, int j) {
            _i = i;
            _j = j;
        }

        inline int getI() const { return _i; }
        inline int getJ() const { return _j; }

    private:
        int _i, _j;
};

class thingy {
        int bogus;
    public:
        thingy(int b) {
            bogus = b;
        }

        int operator[](myindex i);

        int operator()(myindex i, int j);
};

int thingy::operator[](myindex i)
{
    std::cout << "operator[] on index " << i.getI() << ", " << i.getJ() << "\n";
    return bogus;
}

int thingy::operator()(myindex i, int j)
{
    std::cout << "operator() on index " << i.getI() << ", " << i.getJ() << "; " << j << "\n";
    return bogus;
}

int main()
{
    myindex j(1, 5);
    myindex k(3, 4);

    thingy A(14);

    std::cout << A[j] << "\n";
    std::cout << A(k, 5) << "\n";

    return 0;
}
