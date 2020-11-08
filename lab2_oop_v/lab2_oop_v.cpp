// lab2_oop_v.cpp : This file contains the 'main' function. Program execution begins and ends there.

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
using namespace std;
constexpr auto PI = 3.1415926;

class LongNum
{
public:
    LongNum() { this->pos = true; };
    bool operator > (const LongNum& an) const { bool con = (*this <= an); con = !con; return con; };
    bool operator >= (const LongNum& an) const { bool con = (*this < an); con = !con; return con; };
    LongNum(const LongNum& an) { this->fig = an.fig; this->pos = an.pos; };
    const LongNum operator + () const { return LongNum(*this); };
    const LongNum operator - () const { LongNum copy(*this); copy.pos = !copy.pos; return copy; };
    LongNum(const vector<int>& v) { this->fig = v; this->pos = true; };
    bool operator != (const LongNum& an) const { bool con = (*this == an); con = !con; return con; };
    bool operator <= (const LongNum& an) const { bool con1 = (*this < an); bool con2 = (*this == an); return (con1 || con2); };
    LongNum(const string& num)
    {
        string s = num;
        int len = s.length();
        if (len == 0)
        {
            this->pos = true;
        }
        else
        {
            if (s[0] != '-')
            {
                this->pos = true;
            }
            else
            {
                s = s.substr(1);
                this->pos = false;
            }
            int i = len - 1;
            while (i >= 0)
            {
                this->fig.push_back((int)s[i] - '0');
                i -= 1;
            }
            this->r_o();
        }
    };
    bool operator == (const LongNum& an) const {
        if (this->pos == an.pos)
        {
            if (this->fig.empty())
            {
                if (an.fig.empty())
                {
                    return true;
                }
                else if (an.fig.size() == 1)
                {
                    if (an.fig[0] == 0)
                    {
                        return true;
                    }
                    else
                    {
                        return false;
                    }
                }
                else
                {
                    return false;
                }
            }
            if (an.fig.empty())
            {
                if (this->fig.size() == 1)
                {
                    if (this->fig[0] == 0)
                    {
                        return true;
                    }
                    else
                    {
                        return false;
                    }
                    return false;
                }
                else
                {
                    return false;
                }
            }
            if (this->fig.size() != an.fig.size())
            {
                return false;
            }
            for (unsigned int i = 0; i < this->fig.size(); i++)
            {
                if (this->fig[i] != an.fig[i])
                {
                    return false;
                }
            }
            return true;
        }
        else
        {
            return false;
        }
    };
    bool operator < (const LongNum& an) const
    {
        if (*this != an)
        {
            if (!this->pos)
            {
                if (an.pos)
                {
                    return true;
                }
                else
                {
                    bool con = ((-an) < (-(*this)));
                    return con;
                }
            }
            else  if (!an.pos)
            {
                return false;
            }
            else
            {
                if (this->fig.size() == an.fig.size())
                {
                    unsigned int i = this->fig.size() - 1;
                    while (i >= 0)
                    {
                        if (this->fig[i] != an.fig[i])
                        {
                            return (this->fig[i] < an.fig[i]);
                        }
                        i -= 1;
                    }
                    return false;
                }
                else
                {
                    return (this->fig.size() < an.fig.size());
                }
            }
        }
        else
        {
            return false;
        }

    };
    LongNum(const int& n)
    {
        int i = n;
        if (n < 0)
        {
            this->pos = false;
        }
        else
        {
            this->pos = true;
        }
        while (i != 0)
        {
            this->fig.push_back(i % 10);
            i /= 10;
        }
    };
    LongNum& operator = (const LongNum& an) { this->fig = an.fig; this->pos = an.pos; return (*this); };
    LongNum operator + (const LongNum& an) const
    {
        if (!this->pos)
        {
            if (an.pos)
            {
                return (an - (-(*this)));
            }
            else
            {
                return -(-(*this) + (-an));
            }
        }
        else if (!an.pos)
        {
            return ((*this) - (-an));
        }
        vector<int> res;
        int c = 0;
        for (unsigned i = 0; i < max(this->fig.size(), an.fig.size()) || c != 0; ++i) {
            int tot = c;
            tot = tot + (i < this->fig.size() ? this->fig[i] : 0);
            tot = tot + (i < an.fig.size() ? an.fig[i] : 0);
            if (tot <= 9)
            {
                c = 0;
                res.push_back(tot);
            }
            else
            {
                c = 1;
                res.push_back(tot - 10);
            }
        }
        LongNum result(res);
        result.r_o();
        return result;
    };
    LongNum(const long long& n)
    {
        long long i = n;
        if (n < 0)
        {
            this->pos = false;
        }
        else
        {
            this->pos = true;
        }
        while (i != 0)
        {
            this->fig.push_back(i % 10);
            i /= 10;
        }
    };
    LongNum operator - (const LongNum& an) const
    {
        if (!an.pos)
        {
            return (*this + (-an));
        }
        else if (!this->pos)
        {
            return -((-(*this)) + an);
        }
        else if (*this < an)
        {
            return -(an - (*this));
        }
        vector<int> res;
        int c = 0;
        for (unsigned i = 0; i < max(an.fig.size(), this->fig.size()) || c != 0; i++)
        {
            int tot = this->fig[i];
            tot = tot - (c + (i < an.fig.size() ? an.fig[i] : 0));
            if (tot >= 0)
            {
                c = 0;
                res.push_back(tot);
            }
            else
            {
                c = 1;
                res.push_back(tot + 10);
            }
        }
        LongNum r(res);
        r.r_o();
        return r;
    };
    LongNum operator * (const LongNum& an) const
    {
        LongNum res;
        res.fig.resize(this->fig.size() + an.fig.size() + 1, 0);
        for (unsigned i = 0; i < this->fig.size(); i++) 
        {
            int c = 0;
            for (unsigned j = 0; j < an.fig.size() || c != 0; ++j) 
            {
                int tot = res.fig[i + j];
                tot = tot + this->fig[i] * (j < an.fig.size() ? an.fig[j] : 0);
                tot = tot + c;
                res.fig[i + j] = static_cast<int>(tot % 10);
                c = static_cast<int>(tot / 10);
            }
        }
        res.pos = (this->pos == an.pos);
        res.r_o();
        return res;
    };
    void s_r()
    {
        if (this->fig.size() == 0)
        {
            this->fig.push_back(0);
            return;
        }
        this->fig.push_back(this->fig[this->fig.size() - 1]);
        unsigned int i = this->fig.size() - 2;
        while (i > 0)
        {
            this->fig[i] = this->fig[i - 1];
            i -= 1;
        }
        this->fig[0] = 0;
    };
    LongNum operator / (const LongNum& an) const
    {
        LongNum b = an;
        b.pos = true;
        LongNum res, cur;
        res.fig.resize(this->fig.size());
        long long i = static_cast<long long>(this->fig.size()) - 1;
        while (i >= 0)
        {
            cur.s_r();
            cur.fig[0] = this->fig[i];
            cur.r_o();
            int x(0);
            int y(0);
            int z = 10;
            while (y <= z) 
            {
                int m1 = (y + z) / 2;
                LongNum m_(to_string(m1));
                LongNum t = b * m_;
                if (t > cur) 
                {
                    z = m1 - 1;
                }
                else
                {
                    x = m1;
                    y = m1 + 1;
                }
            }
            res.fig[i] = x;
            LongNum v(to_string(x));
            cur = cur - b * z;
            i -= 1;
        }
        res.pos = (this->pos == an.pos);
        res.r_o();
        return res;
    };
    LongNum operator % (const LongNum& mod) const
    {
        if (*this < mod)
        {
            return *this;
        }
        else if (*this == mod)
        {
            LongNum nul("0");
            return nul;
        }
        LongNum res;
        LongNum r = *this / mod;
        res = *this - r * mod;
        return res;
    };
    void s_l() {
        if (this->fig.size() == 0)
        {
            this->fig.push_back(0);
            return;
        }

        this->fig.erase(fig.begin());
    };
    LongNum gcd(const LongNum& an) {
        LongNum a = *this, b = an;
        while (a != b) 
        {
            if (a > b) 
            {
                LongNum f = a;
                a = b;
                b = f;
            }
            b = b - a;
        }
        return a;
    }
    LongNum random() {
        LongNum one("1");
        LongNum substract = *this - one;
        vector<int> u = substract.fig;
        int a;
        for (int i = 0; i < u.size(); i++) 
        {
            if (u[i] == 0)
            {
                continue;
            }
            a = 1 + rand() % u[i];
            u[i] = a;
        }
        LongNum res(u);
        res.r_o();
        return res;
    }
    vector<LongNum> fact()
    {
        LongNum x = (*this);
        vector<LongNum> res;
        LongNum n = x / 2;
        LongNum one("1");
        LongNum nul("0");
        LongNum i("2");
        for (i; i <= n; i = i + one) 
        {
            while (x % i == nul) 
            {
                res.push_back(i);
                x = x / i;
            }
        }
        return res;
    };
    void out()
    {
        if (!this->pos)
        {
            cout << "-";
        }
        int i = (int)this->fig.size() - 1;
        while (i >= 0)
        {
            cout << this->fig[i];
            i -= 1;
        }
    }
    void r_o()
    {
        while (this->fig.size() > 1 && this->fig.back() == 0) 
        {
            this->fig.pop_back();
        }
        if (this->fig.size() == 1 && this->fig[0] == 0)
        {
            this->pos = true;
        }
    };
    bool od() const
    {
        if (this->fig.size() == 0)
        {
            return false;
        }
        bool r = (this->fig[0] & 1);
        return r;
    }
    LongNum pow(LongNum n) const 
    {
        LongNum a(*this);
        LongNum res(1);
        while (n != 0) 
        {
            if (n.od())
            {
                res = res * a;
            }
            a = a * a;
            n = n / LongNum(2);
        }
        return res;
    }
    LongNum sq_bin()
    {
        LongNum one("1");
        LongNum two("2");
        LongNum o("0");
        LongNum t = (*this) / two;
        while (o <= t) 
        {
            LongNum m = (t + o) / two;
            LongNum y = m * m;
            if (y <= (*this)) 
            {
                o = m + one;
            }
            else
            {
                t = m - one;
            }
        }
        return o - one;
    };
    LongNum to_bin() const
    {
        LongNum z("0");
        vector<int> bin_repr;
        LongNum con = *this;

        while (con > z) 
        {
            bin_repr.push_back(con.fig[0] % 2);
            con = con / 2;
        }
        if (!this->pos) 
        {
            bin_repr.at(0) = 0;
        }
        return LongNum(bin_repr);
    };
    vector<int> sum(const vector<int>& x, const vector<int>& y)
    {
        int c = 0;
        vector<int> res;
        for (unsigned int i = 0; i < max(x.size(), y.size()) || c != 0; ++i) 
        {
            int tot = c;
            tot = tot + (i < x.size() ? x[i] : 0);
            tot = tot + (i < y.size() ? y[i] : 0);
            if (tot > 9) 
            {
                c = 1;
                res.push_back(tot - 10);
            }
            else 
            {
                c = 0;
                res.push_back(tot);
            }
        }
        while (res.size() > 1 && res.back() == 0)
        {
            res.pop_back();
        }
        return res;
    };
    vector<int> substract(const vector<int>& x, const vector<int>& y)
    {
        int c = 0;
        vector<int> res;
        for (unsigned i = 0; i < max(x.size(), y.size()) || c != 0; i++) 
        {
            int tot = x[i];
            tot = tot - (c + (i < y.size() ? y[i] : 0));
            if (tot >= 0) 
            {
                c = 0;
                res.push_back(tot);
            }
            else 
            {
                c = 1;
                res.push_back(tot + 10);
            }
        }
        while (res.size() > 1 && res.back() == 0)
        {
            res.pop_back();
        }
        return res;
    };
    vector<int> multiplication(const vector<int>& x, const vector<int>& y)
    {
        vector<int> res(x.size() + y.size() + 1, 0);
        for (unsigned i = 0; i < x.size(); i++) 
        {
            int c = 0;
            for (unsigned j = 0; j < y.size() || c != 0; ++j) 
            {
                int cur = res[i + j];
                cur = cur + x[i] * (j < y.size() ? y[j] : 0);
                cur = cur + c;
                res[i + j] = static_cast<int>(cur % 10);
                c = static_cast<int>(cur / 10);
            }
        }
        while (res.size() > 1 && res.back() == 0)
        {
            res.pop_back();
        }
        return res;
    };
    vector<int> pow(const vector<int>& x, long long n) {
        vector<int> y = x;
        vector<int> res = {1};
        while (n != 0) 
        {
            if (n % 2 == 1)
            {
                res = multiplication(res, y);
            }
            y = multiplication(y, y);
            n /= 2;
        }
        return res;
    };
    void ft(vector<int>& x, bool in)
    {
        int n = (int)x.size();
        if (n == 1)
        {
            return;
        }
        vector<int> a0(n / 2);
        vector<int> a1(n / 2);
        for (int i = 0, j = 0; i < n; i += 2, j++) 
        {
            a0[j] = x[i];
            a1[j] = x[i + 1];
        }
        ft(a0, in);
        ft(a1, in);
        double ang = 2 * PI / n * (in ? -1 : 1);
        int e(1);
        for (int i = 0; i < n / 2; i++) 
        {
            x[i] = a0[i] + e * a1[i];
            x[i + n / 2] = a0[i] - e * a1[i];
            if (in)
            {
                x[i] = x[i] / 2;
                x[i + n / 2] = x[i + n / 2] / 2;
            }
            e *= e;
        }
    };
    LongNum karatsuba(const LongNum& x)
    {
        vector<int> a = x.fig;
        vector<int> b = this->fig;
        int len = max(a.size(), b.size());
        int len1 = len;
        int len2 = len;
        while (len1 & (len1 - 1))
        {
            len1 += 1;
        }
        while (len2 & (len2 - 1))
        {
            len2 += 1;
        }
        b.resize(len2);
        vector<int> res = v_k(a, b);
        LongNum res_(res);
        return res_;
    };
    vector<int> v_k(const vector<int>& a, const vector<int>& b)
    {
        int len = a.size();
        if (len <= 216) 
        {
            return multiplication(a, b);
        }
        long long k = len / 2;
        string dec = "10";
        vector<int> aq{a.begin() + k, a.end()};
        vector<int> bq{b.begin() + k, b.end()};
        vector<int> as{a.begin(), a.begin() + k};
        vector<int> bs{b.begin(), b.begin() + k};
        vector<int> c1 = sum(as, aq);
        vector<int> c2 = sum(bs, bq);
        vector<int> d1 = v_k(aq, bq);
        vector<int> d2 = v_k(as, bs);
        vector<int> d3 = v_k(c1, c2);
        vector<int> res;
        int i = dec.size() - 1;
        while (i >= 0)
        {
            res.push_back(dec[i] - 48);
            i -= 1;
        }
        vector<int> dec1 = res;
        vector<int> e1 = pow(dec1, len);
        vector<int> e2 = pow(dec1, len / 2);
        vector<int> res1 = multiplication(d1, e1);
        vector<int> res2 = substract(d3, d1);
        res2 = substract(res2, d2);
        res2 = multiplication(res2, e2);
        vector<int> res_ = sum(res1, res2);
        res_ = sum(res_, d2);
        return res_;
    }
    LongNum shonhage(const LongNum& num)
    {
        int len = max(this->fig.size(), num.fig.size());
        if (len <= 216) 
        {
            return *this * num;
        }
        unsigned int p0 = 1;
        LongNum n("2");
        while (*this > n || num > n) 
        {
            p0++;
            n = n + n;
        }
        int j = ceil(double(p0 - 8) / 18);
        LongNum n1((int)round(std::pow(2, 6 * j - 1) - 1));
        LongNum n2((int)round(std::pow(2, 6 * j + 1) - 1));
        LongNum n3((int)round(std::pow(2, 6 * j + 2) - 1));
        LongNum n4((int)round(std::pow(2, 6 * j + 3) - 1));
        LongNum n5((int)round(std::pow(2, 6 * j + 5) - 1));
        LongNum n6((int)round(std::pow(2, 6 * j + 7) - 1));
        LongNum k1 = *this % n1;
        LongNum k2 = *this % n2;
        LongNum k3 = *this % n3;
        LongNum k4 = *this % n4;
        LongNum k5 = *this % n5;
        LongNum k6 = *this % n6;
        LongNum p1 = num % n1;
        LongNum p2 = num % n2;
        LongNum p3 = num % n3;
        LongNum p4 = num % n4;
        LongNum p5 = num % n5;
        LongNum p6 = num % n6;
        LongNum o1 = k1.karatsuba(p1) % n1;
        LongNum o2 = k2.karatsuba(p2) % n2;
        LongNum o3 = k3.karatsuba(p3) % n3;
        LongNum o4 = k4.karatsuba(p4) % n4;
        LongNum o5 = k5.karatsuba(p5) % n5;
        LongNum o6 = k6.karatsuba(p6) % n6;
        LongNum o = (num > * this ? num : *this);
        o = o - (o % n6) + o6;
        while (o % n1 != o1 || o % n2 != o2 || o % n3 != o3 || o % n4 != o4 || o % n5 != o || o % n6 != o6) 
        {
            o = o + n6;
        }
        return o;
    };
    LongNum inv_num() const
    {
        LongNum zero("0");
        LongNum one("1");
        LongNum two("2");
        LongNum four("4");
        if (*this == zero) 
        { 
            return zero; 
        }
        LongNum q((-1) + (int)this->to_bin().fig.size());
        LongNum f = two.pow(q);
        LongNum s = f;
    it:
        f = two * f - (*this) * ((f * f) / two.pow(q)) / two.pow(q);
        if (f <= s)
        {
            goto res;
        }
        s = f;
        goto it;
    res:
        LongNum y(four.pow(q) - (*this) * f);
        while (y < zero) 
        {
            f = f - one;
            y = y + *this;
        }
        return f;
    };
    LongNum div_by(const LongNum& num)
    {
        LongNum one("1");
        long long shi = 2 * (static_cast<long long> ((num - one).to_bin().fig.size()));
        LongNum inverse_divider = num.inv_num();
        LongNum result = inverse_divider.karatsuba(*this);
        for (long long i = 1; i <= shi; i *= 2)
        {
            result.s_l();
        }
        return *this / num;
    };
    bool sol_stras_prim_test()
    {
        LongNum two("2");
        LongNum one("1");
        LongNum t = *this;
        const int m = 5;
        if (*this <= two) 
        { 
            return true; 
        }
        if (*this % two == 0)
        { 
            return false; 
        }
        vector<LongNum> v = t.fact();
        LongNum x = t - one;
        for (int i = 1; i <= m; i++) 
        {
            LongNum a = t.random();
            LongNum g = a.gcd(t);
            if (g > one) 
            { 
                return false; 
            }
            LongNum c = a.pow(x);
            c = c % t;
            LongNum r = jac_sym(a, t, v);
            r = r % t;
            if (c != r) 
            { 
                return false; 
            }
        }
        return true;
    };
    LongNum jac_sym(LongNum a, LongNum b, const vector<LongNum>& b_fact)
    {
        LongNum res("1");
        for (long long i = 0; i < b_fact.size(); i++)
        {
            LongNum x = leg_sym(a, b_fact[i]);
            res = res * x;
        }
        return res;
    };
    LongNum leg_sym(const LongNum& a, const LongNum& b)
    {
        LongNum check = a % b;
        LongNum nul("0");
        LongNum one("1");
        LongNum two("2");
        if (check == nul) 
        {
            return nul;
        }
        LongNum p = a % b;
        for (LongNum i("1"); i <= b; i = i + one) 
        {
            LongNum z = i * i;
            z = z % b;
            if (z == p) 
            {
                return one;
            }
        }
        LongNum neq("-1");
        return neq;
    };
    bool mil_rab_p_test()
    {
        LongNum two("2");
        LongNum one("1");
        LongNum z = (*this) - one;
        LongNum s("1");
        LongNum t = z;
        if (*this <= two) 
        { 
            return true; 
        }
        if (*this % two == 0) 
        { 
            return false; 
        }
        while (true) 
        {
            t = t / two;
            if (t % two == 1) { break; }
            s = s + one;
        }
        const int k = 5;
        for (int i = 1; i <= k; i++) 
        {
            LongNum b = z.random();
            LongNum o = b.pow(t);
            o = o % (*this);
            if (o == one || o == z)
            { 
                continue; 
            }
            for (LongNum i("1"); i < s; i = i + one) 
            {
                o = o * o;
                o = o % (*this);
                if (o == one) 
                { 
                    return false; 
                }
                else if (o == z)
                {
                    break;
                }
            }
            if (o == z) 
            { 
                continue; 
            }
            return false;
        }
        return true;
    };
private:
    vector<int> fig;
    bool pos;
};

int main()
{
    string x;
    string y;
    cout << "x = ";
    cin >> x;
    cout << "y = ";
    cin >> y;
    LongNum num1(x);
    LongNum num2(y);
    LongNum res = num1.karatsuba(num2);
    cout << "Karatsuba - ";
    res.out();
    cout << endl;
    cout << "Strasen-Shonhage mult - ";
    res = num1.shonhage(num2);
    res.out();
    cout << endl;
    cout << "Inv number - ";
    res = res.inv_num();
    res.out();
    cout << endl;
    cout << "Integer division - ";
    res = num1.div_by(num2);
    res.out();
    cout << endl;
    cout << "Solovay-Strasen test  - ";
    if (num1.sol_stras_prim_test())
    {
        cout << "Prime" << endl;
    }
    else
    {
        cout << "Not prime" << endl;
    }
    cout << "Miller-Rabin test - ";
    if (num1.mil_rab_p_test())
    {
        cout << "Prime" << endl;
    }
    else
    {
        cout << "Not prime" << endl;
    }
    return 0;
}