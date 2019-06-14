#ifndef UNITS_GCT
#define UNITS_GCT

namespace units{
    constexpr double mtopc = 3.24078e-17;
    constexpr double pctom = 3.08568e16;

    constexpr double MeVc2toKg = 1.78266e-30;
    constexpr double KgtoMeVc2 = 5.60959e29;

    constexpr double Ctoe = 6.24151e18;
    constexpr double etoC = 1.60218e-19;

    constexpr double TtomG = 1e10;
    constexpr double mGtoT = 1e-10;

    constexpr double Ytos = 31557600.0;
    constexpr double stoY = 3.17098e-8;
}

/********** DISTANCE **********/

class Meter{
    public:

        double value;

      // Conversion operator to allow [class] a = Meter b
        operator double() {return this->value;}
        operator Parsec() {return this->value * units::mtopc;}
      
      // Assignement operator to allow Meter a = [class] b
        void operator=(const double &rhs)   {this->value = rhs;}
        void operator=(const int &rhs)      {this->value = static_cast<double>(rhs);}
        void operator=(const Parsec &rhs)   {this->value = units::pctom * rhs.value;}

      // MDAS

        Meter operator+(const Meter &rhs) const    { return { value + rhs.value }; }
        Meter operator-(const Meter &rhs) const    { return { value - rhs.value }; }
        Meter operator*(const Meter &rhs) const    { return { value * rhs.value }; }
        Meter operator/(const Meter &rhs) const    { return { value / rhs.value }; }

        Meter operator+(const int &rhs)    { return { value + static_cast<double>(rhs) }; }
        Meter operator+(const double &rhs) { return { value + rhs }; }
        Meter operator+(const Parsec &rhs) { return { value + rhs.value * units::pctom}; }

        Meter operator-(const int &rhs)    { return { value - static_cast<double>(rhs) }; }
        Meter operator-(const double &rhs) { return { value - rhs }; }
        Meter operator-(const Parsec &rhs) { return { value - rhs.value * units::pctom}; }

        Meter operator*(const int &rhs)    { return { value * static_cast<double>(rhs) }; }
        Meter operator*(const double &rhs) { return { value * rhs }; }
        Meter operator*(const Parsec &rhs) { return { value * rhs.value * units::pctom}; }

        Meter operator/(const int &rhs)    { return { value / static_cast<double>(rhs) }; }
        Meter operator/(const double &rhs) { return { value / rhs }; }
        Meter operator/(const Parsec &rhs) { return { value / (rhs.value * units::pctom)}; }

      // Constructors from different classes. Allows std::array<Meter, n> a = { 1, 2, 3, ... }
        Meter()                     {this->value = 0;}
        Meter(const double &doub)   {this->value = doub;}
        Meter(const int &iint)      {this->value = static_cast<double>(iint);}
        Meter(const Parsec &pc)     {this->value = units::pctom * pc.value;}
        Meter(const Meter &m)       {this->value = m.value;}
};


class Parsec{
    public:

        double value;

      // Conversion operator to allow [class] a = Parsec b
        operator double()   {return this->value;}
        operator Meter()    {return this->value * units::pctom;}

      // Assignement operator to allow Parsec a = [class] b
        void operator=(const double &rhs)   {this->value = rhs;}
        void operator=(const int &rhs)      {this->value = static_cast<double>(rhs);}
        void operator=(const Meter &rhs)    {this->value = units::mtopc * rhs.value;}

      // MDAS
        Parsec operator+(const Parsec &rhs) const    { return { value + rhs.value }; }
        Parsec operator-(const Parsec &rhs) const    { return { value - rhs.value }; }
        Parsec operator*(const Parsec &rhs) const    { return { value * rhs.value }; }
        Parsec operator/(const Parsec &rhs) const    { return { value / rhs.value }; }

        Parsec operator+(const int &rhs)    { return { value + static_cast<double>(rhs) }; }
        Parsec operator+(const double &rhs) { return { value + rhs }; }
        Parsec operator+(const Meter &rhs)  { return { value + rhs.value * units::mtopc}; }

        Parsec operator-(const int &rhs)    { return { value - static_cast<double>(rhs) }; }
        Parsec operator-(const double &rhs) { return { value - rhs }; }
        Parsec operator-(const Meter &rhs)  { return { value - rhs.value * units::mtopc}; }

        Parsec operator*(const int &rhs)    { return { value * static_cast<double>(rhs) }; }
        Parsec operator*(const double &rhs) { return { value * rhs }; }
        Parsec operator*(const Meter &rhs)  { return { value * rhs.value * units::mtopc}; }

        Parsec operator/(const int &rhs)    { return { value / static_cast<double>(rhs) }; }
        Parsec operator/(const double &rhs) { return { value / rhs }; }
        Parsec operator/(const Meter &rhs)  { return { value / (rhs.value * units::mtopc)}; }

      // Constructors from different classes. Allows std::array<Parsec, n> a = { 1, 2, 3, ... }
        Parsec()                     {this->value = 0;}
        Parsec(const double &doub)   {this->value = doub;}
        Parsec(const int &iint)      {this->value = static_cast<double>(iint);}
        Parsec(const Meter &m)       {this->value = units::mtopc * m.value;}
        Parsec(const Parsec &pc)     {this->value = pc.value;}
};


/********** MASS **********/

class MeVc2{
    public:

        double value;

      // Conversion operator to allow [class] a = MeVc2 b
        operator double()   {return this->value;}
        operator Kg()       {return this->value * units::MeVc2toKg;}
      
      // Assignement operator to allow MeVc2 a = [class] b
        void operator=(const double &rhs)   {this->value = rhs;}
        void operator=(const int &rhs)      {this->value = static_cast<double>(rhs);}
        void operator=(const Kg &rhs)       {this->value = units::KgtoMeVc2 * rhs.value;}

      // MDAS
        MeVc2 operator+(const MeVc2 &rhs) const    { return { value + rhs.value }; }
        MeVc2 operator-(const MeVc2 &rhs) const    { return { value - rhs.value }; }
        MeVc2 operator*(const MeVc2 &rhs) const    { return { value * rhs.value }; }
        MeVc2 operator/(const MeVc2 &rhs) const    { return { value / rhs.value }; }

        MeVc2 operator+(const int &rhs)    { return { value + static_cast<double>(rhs) }; }
        MeVc2 operator+(const double &rhs) { return { value + rhs }; }
        MeVc2 operator+(const Kg &rhs)     { return { value + rhs.value * units::KgtoMeVc2}; }

        MeVc2 operator-(const int &rhs)    { return { value - static_cast<double>(rhs) }; }
        MeVc2 operator-(const double &rhs) { return { value - rhs }; }
        MeVc2 operator-(const Kg &rhs)     { return { value - rhs.value * units::KgtoMeVc2}; }

        MeVc2 operator*(const int &rhs)    { return { value * static_cast<double>(rhs) }; }
        MeVc2 operator*(const double &rhs) { return { value * rhs }; }
        MeVc2 operator*(const Kg &rhs)     { return { value * rhs.value * units::KgtoMeVc2}; }

        MeVc2 operator/(const int &rhs)    { return { value / static_cast<double>(rhs) }; }
        MeVc2 operator/(const double &rhs) { return { value / rhs }; }
        MeVc2 operator/(const Kg &rhs)     { return { value / (rhs.value * units::KgtoMeVc2)}; }

      // Constructors from different classes. Allows std::array<MeVc2, n> a = { 1, 2, 3, ... }
        MeVc2()                     {this->value = 0;}
        MeVc2(const double &doub)   {this->value = doub;}
        MeVc2(const int &iint)      {this->value = static_cast<double>(iint);}
        MeVc2(const Kg &kg)         {this->value = units::KgtoMeVc2 * kg.value;}
        MeVc2(const MeVc2 &m)       {this->value = m.value;}
};

class Kg{
    public:

        double value;

      // Conversion operator to allow [class] a = Kg b
        operator double()   {return this->value;}
        operator Kg()       {return this->value * units::MeVc2toKg;}

      // Assignement operator to allow MeVc2 a = [class] b
        void operator=(const double &rhs)   {this->value = rhs;}
        void operator=(const int &rhs)      {this->value = static_cast<double>(rhs);}
        void operator=(const MeVc2 &rhs)    {this->value = units::MeVc2toKg * rhs.value;}

      // MDAS
        Kg operator+(const Kg &rhs) const    { return { value + rhs.value }; }
        Kg operator-(const Kg &rhs) const    { return { value - rhs.value }; }
        Kg operator*(const Kg &rhs) const    { return { value * rhs.value }; }
        Kg operator/(const Kg &rhs) const    { return { value / rhs.value }; }
        
        Kg operator+(const int &rhs)    { return { value + static_cast<double>(rhs) }; }
        Kg operator+(const double &rhs) { return { value + rhs }; }
        Kg operator+(const MeVc2 &rhs)  { return { value + rhs.value * units::MeVc2toKg}; }

        Kg operator-(const int &rhs)    { return { value - static_cast<double>(rhs) }; }
        Kg operator-(const double &rhs) { return { value - rhs }; }
        Kg operator-(const MeVc2 &rhs)  { return { value - rhs.value * units::MeVc2toKg}; }

        Kg operator*(const int &rhs)    { return { value * static_cast<double>(rhs) }; }
        Kg operator*(const double &rhs) { return { value * rhs }; }
        Kg operator*(const MeVc2 &rhs)  { return { value * rhs.value * units::MeVc2toKg}; }

        Kg operator/(const int &rhs)    { return { value / static_cast<double>(rhs) }; }
        Kg operator/(const double &rhs) { return { value / rhs }; }
        Kg operator/(const MeVc2 &rhs)  { return { value / (rhs.value * units::MeVc2toKg)}; }

      // Extra MDAS
        double operator/(const Coulomb &q) const { return value / q.value; }

      // Constructors from different classes. Allows std::array<Kg, n> a = { 1, 2, 3, ... }
        Kg()                     {this->value = 0;}
        Kg(const double &doub)   {this->value = doub;}
        Kg(const int &iint)      {this->value = static_cast<double>(iint);}
        Kg(const MeVc2 &m)       {this->value = units::MeVc2toKg * m.value;}
        Kg(const Kg &kg)         {this->value = kg.value;}
};


/********** CHARGE **********/

class e_charge{
    public:

        double value;

      // Conversion operator to allow [class] a = e_charge b
        operator double()   {return this->value;}
        operator Coulomb()  {return this->value * units::etoC;}

      // Assignement operator to allow e_charge a = [class] b
        void operator=(const double &rhs)   {this->value = rhs;}
        void operator=(const int &rhs)      {this->value = static_cast<double>(rhs);}
        void operator=(const Coulomb &rhs)  {this->value = units::Ctoe * rhs.value;}

      // MDAS
        e_charge operator+(const e_charge &rhs) const    { return { value + rhs.value }; }
        e_charge operator-(const e_charge &rhs) const    { return { value - rhs.value }; }
        e_charge operator*(const e_charge &rhs) const    { return { value * rhs.value }; }
        e_charge operator/(const e_charge &rhs) const    { return { value / rhs.value }; }

        e_charge operator+(const int &rhs)      { return { value + static_cast<double>(rhs) }; }
        e_charge operator+(const double &rhs)   { return { value + rhs }; }
        e_charge operator+(const Coulomb &rhs)  { return { value + rhs.value * units::Ctoe}; }

        e_charge operator-(const int &rhs)      { return { value - static_cast<double>(rhs) }; }
        e_charge operator-(const double &rhs)   { return { value - rhs }; }
        e_charge operator-(const Coulomb &rhs)  { return { value - rhs.value * units::Ctoe}; }

        e_charge operator*(const int &rhs)      { return { value * static_cast<double>(rhs) }; }
        e_charge operator*(const double &rhs)   { return { value * rhs }; }
        e_charge operator*(const Coulomb &rhs)  { return { value * rhs.value * units::Ctoe}; }

        e_charge operator/(const int &rhs)      { return { value / static_cast<double>(rhs) }; }
        e_charge operator/(const double &rhs)   { return { value / rhs }; }
        e_charge operator/(const Coulomb &rhs)  { return { value / (rhs.value * units::Ctoe)}; }

      // Constructors from different classes. Allows std::array<e_charge, n> a = { 1, 2, 3, ... }
        e_charge()                     {this->value = 0;}
        e_charge(const double &doub)   {this->value = doub;}
        e_charge(const int &iint)      {this->value = static_cast<double>(iint);}
        e_charge(const Coulomb &c)     {this->value = units::Ctoe * c.value;}
        e_charge(const e_charge &e)    {this->value = e.value;}

};

class Coulomb{
    public:

        double value;

      // Conversion operator to allow [class] a = Coulomb b
        operator double()   {return this->value;}
        operator e_charge() {return this->value * units::Ctoe;}

      // Assignement operator to allow Coulomb a = [class] b
        void operator=(const double &rhs)   {this->value = rhs;}
        void operator=(const int &rhs)      {this->value = static_cast<double>(rhs);}
        void operator=(const Coulomb &rhs)  {this->value = units::etoC * rhs.value;}

      // MDAS
        Coulomb operator+(const Coulomb &rhs) const    { return { value + rhs.value }; }
        Coulomb operator-(const Coulomb &rhs) const    { return { value - rhs.value }; }
        Coulomb operator*(const Coulomb &rhs) const    { return { value * rhs.value }; }
        Coulomb operator/(const Coulomb &rhs) const    { return { value / rhs.value }; }

        Coulomb operator+(const int &rhs)       { return { value + static_cast<double>(rhs) }; }
        Coulomb operator+(const double &rhs)    { return { value + rhs }; }
        Coulomb operator+(const e_charge &rhs)  { return { value + rhs.value * units::etoC}; }

        Coulomb operator-(const int &rhs)       { return { value - static_cast<double>(rhs) }; }
        Coulomb operator-(const double &rhs)    { return { value - rhs }; }
        Coulomb operator-(const e_charge &rhs)  { return { value - rhs.value * units::etoC}; }

        Coulomb operator*(const int &rhs)       { return { value * static_cast<double>(rhs) }; }
        Coulomb operator*(const double &rhs)    { return { value * rhs }; }
        Coulomb operator*(const e_charge &rhs)  { return { value * rhs.value * units::etoC}; }

        Coulomb operator/(const int &rhs)       { return { value / static_cast<double>(rhs) }; }
        Coulomb operator/(const double &rhs)    { return { value / rhs }; }
        Coulomb operator/(const e_charge &rhs)  { return { value / (rhs.value * units::etoC)}; }

      // Extra MDAS
        friend double operator/(const double &lhs, const Coulomb &rhs);

      // Constructors from different classes. Allows std::array<Coulomb, n> a = { 1, 2, 3, ... }
        Coulomb()                     {this->value = 0;}
        Coulomb(const double &doub)   {this->value = doub;}
        Coulomb(const int &iint)      {this->value = static_cast<double>(iint);}
        Coulomb(const e_charge &e)    {this->value = units::etoC * e.value;}
        Coulomb(const Coulomb &c)     {this->value = c.value;}
};

        double operator/(const double &lhs, const Coulomb &rhs){
            return lhs/rhs.value;
        }

/********** MAGNETIC FIELD **********/

class microGauss{
    public:

        double value;

      // Conversion operator to allow [class] a = microGauss b
        operator double()   {return this->value;}
        operator Tesla()    {return this->value * units::mGtoT;}

      // Assignement operator to allow microGauss a = [class] b
        void operator=(const double &rhs)   {this->value = rhs;}
        void operator=(const int &rhs)      {this->value = static_cast<double>(rhs);}
        void operator=(const Tesla &rhs)    {this->value = units::TtomG * rhs.value;}

      // MDAS
        microGauss operator+(const microGauss &rhs) const    { return { value + rhs.value }; }
        microGauss operator-(const microGauss &rhs) const    { return { value - rhs.value }; }
        microGauss operator*(const microGauss &rhs) const    { return { value * rhs.value }; }
        microGauss operator/(const microGauss &rhs) const    { return { value / rhs.value }; }

        microGauss operator+(const int &rhs)    { return { value + static_cast<double>(rhs) }; }
        microGauss operator+(const double &rhs) { return { value + rhs }; }
        microGauss operator+(const Tesla &rhs)  { return { value + rhs.value * units::TtomG}; }

        microGauss operator-(const int &rhs)    { return { value - static_cast<double>(rhs) }; }
        microGauss operator-(const double &rhs) { return { value - rhs }; }
        microGauss operator-(const Tesla &rhs)  { return { value - rhs.value * units::TtomG}; }

        microGauss operator*(const int &rhs)    { return { value * static_cast<double>(rhs) }; }
        microGauss operator*(const double &rhs) { return { value * rhs }; }
        microGauss operator*(const Tesla &rhs)  { return { value * rhs.value * units::TtomG}; }

        microGauss operator/(const int &rhs)    { return { value / static_cast<double>(rhs) }; }
        microGauss operator/(const double &rhs) { return { value / rhs }; }
        microGauss operator/(const Tesla &rhs)  { return { value / (rhs.value * units::TtomG)}; }

      // Constructors from different classes. Allows std::array<microGauss, n> a = { 1, 2, 3, ... }
        microGauss()                     {this->value = 0;}
        microGauss(const double &doub)   {this->value = doub;}
        microGauss(const int &iint)      {this->value = static_cast<double>(iint);}
        microGauss(const Tesla &T)       {this->value = units::TtomG * T.value;}
        microGauss(const microGauss &mG) {this->value = mG.value;}
};

class Tesla{
    public:

        double value;

      // Conversion operator to allow [class] a = Tesla b
        operator double()       {return this->value;}
        operator microGauss()   {return this->value * units::TtomG;}

      // Assignement operator to allow Tesla a = [class] b
        void operator=(const double &rhs)       {this->value = rhs;}
        void operator=(const int &rhs)          {this->value = static_cast<double>(rhs);}
        void operator=(const microGauss &rhs)   {this->value = units::mGtoT * rhs.value;}

      // MDAS
        Tesla operator+(const Tesla &rhs) const    { return { value + rhs.value }; }
        Tesla operator-(const Tesla &rhs) const    { return { value - rhs.value }; }
        Tesla operator*(const Tesla &rhs) const    { return { value * rhs.value }; }
        Tesla operator/(const Tesla &rhs) const    { return { value / rhs.value }; }

        Tesla operator+(const int &rhs)    { return { value + static_cast<double>(rhs) }; }
        Tesla operator+(const double &rhs) { return { value + rhs }; }
        Tesla operator+(const microGauss &rhs)  { return { value + rhs.value * units::mGtoT}; }

        Tesla operator-(const int &rhs)    { return { value - static_cast<double>(rhs) }; }
        Tesla operator-(const double &rhs) { return { value - rhs }; }
        Tesla operator-(const microGauss &rhs)  { return { value - rhs.value * units::mGtoT}; }

        Tesla operator*(const int &rhs)    { return { value * static_cast<double>(rhs) }; }
        Tesla operator*(const double &rhs) { return { value * rhs }; }
        Tesla operator*(const microGauss &rhs)  { return { value * rhs.value * units::mGtoT}; }

        Tesla operator/(const int &rhs)    { return { value / static_cast<double>(rhs) }; }
        Tesla operator/(const double &rhs) { return { value / rhs }; }
        Tesla operator/(const microGauss &rhs)  { return { value / (rhs.value * units::mGtoT)}; }

      // Constructors from different classes. Allows std::array<Tesla, n> a = { 1, 2, 3, ... }
        Tesla()                     {this->value = 0;}
        Tesla(const double &doub)   {this->value = doub;}
        Tesla(const int &iint)      {this->value = static_cast<double>(iint);}
        Tesla(const microGauss &mG) {this->value = units::mGtoT * mG.value;}
        Tesla(const Tesla &T)       {this->value = T.value;}
};


/********** ELECTRIC FIELD **********/

class Vpm{
    public:

        double value;

      // Conversion operator to allow [class] a = Vpm b
        operator double()       {return this->value;}

      // Assignement operator to allow Vpm a = [class] b
        void operator=(const double &rhs)       {this->value = rhs;}
        void operator=(const int &rhs)          {this->value = static_cast<double>(rhs);}

      // MDAS
        Vpm operator+(const Vpm &rhs) const    { return { value + rhs.value }; }
        Vpm operator-(const Vpm &rhs) const    { return { value - rhs.value }; }
        Vpm operator*(const Vpm &rhs) const    { return { value * rhs.value }; }
        Vpm operator/(const Vpm &rhs) const    { return { value / rhs.value }; }

        Vpm operator+(const int &rhs)    { return { value + static_cast<double>(rhs) }; }
        Vpm operator+(const double &rhs) { return { value + rhs }; }

        Vpm operator-(const int &rhs)    { return { value - static_cast<double>(rhs) }; }
        Vpm operator-(const double &rhs) { return { value - rhs }; }

        Vpm operator*(const int &rhs)    { return { value * static_cast<double>(rhs) }; }
        Vpm operator*(const double &rhs) { return { value * rhs }; }

        Vpm operator/(const int &rhs)    { return { value / static_cast<double>(rhs) }; }
        Vpm operator/(const double &rhs) { return { value / rhs }; }
    
      // Constructors from different classes. Allows std::array<Tesla, n> a = { 1, 2, 3, ... }
        Vpm()                     {this->value = 0;}
        Vpm(const double &doub)   {this->value = doub;}
        Vpm(const int &iint)      {this->value = static_cast<double>(iint);}
};

/********** VELOCITY **********/

class mps{
    public:

        double value;

      // Conversion operator to allow [class] a = mps b
        operator double()       {return this->value;}

      // Assignement operator to allow mps a = [class] b
        void operator=(const double &rhs)       {this->value = rhs;}
        void operator=(const int &rhs)          {this->value = static_cast<double>(rhs);}

      // MDAS
        mps operator+(const mps &rhs) const    { return { value + rhs.value }; }
        mps operator-(const mps &rhs) const    { return { value - rhs.value }; }
        mps operator*(const mps &rhs) const    { return { value * rhs.value }; }
        mps operator/(const mps &rhs) const    { return { value / rhs.value }; }

        mps operator+(const int &rhs)    { return { value + static_cast<double>(rhs) }; }
        mps operator+(const double &rhs) { return { value + rhs }; }

        mps operator-(const int &rhs)    { return { value - static_cast<double>(rhs) }; }
        mps operator-(const double &rhs) { return { value - rhs }; }

        mps operator*(const int &rhs)    { return { value * static_cast<double>(rhs) }; }
        mps operator*(const double &rhs) { return { value * rhs }; }

        mps operator/(const int &rhs)    { return { value / static_cast<double>(rhs) }; }
        mps operator/(const double &rhs) { return { value / rhs }; }
    
      // Constructors from different classes. Allows std::array<Tesla, n> a = { 1, 2, 3, ... }
        mps()                     {this->value = 0;}
        mps(const double &doub)   {this->value = doub;}
        mps(const int &iint)      {this->value = static_cast<double>(iint);}
};

/********** ACCELERATION **********/

class mp_ss{
    public:

        double value;

      // Conversion operator to allow [class] a = mps b
        operator double()       {return this->value;}

      // Assignement operator to allow mps a = [class] b
        void operator=(const double &rhs)       {this->value = rhs;}
        void operator=(const int &rhs)          {this->value = static_cast<double>(rhs);}

      // MDAS
        mp_ss operator+(const mp_ss &rhs) const    { return { value + rhs.value }; }
        mp_ss operator-(const mp_ss &rhs) const    { return { value - rhs.value }; }
        mp_ss operator*(const mp_ss &rhs) const    { return { value * rhs.value }; }
        mp_ss operator/(const mp_ss &rhs) const    { return { value / rhs.value }; }

        mp_ss operator+(const int &rhs)    { return { value + static_cast<double>(rhs) }; }
        mp_ss operator+(const double &rhs) { return { value + rhs }; }

        mp_ss operator-(const int &rhs)    { return { value - static_cast<double>(rhs) }; }
        mp_ss operator-(const double &rhs) { return { value - rhs }; }

        mp_ss operator*(const int &rhs)    { return { value * static_cast<double>(rhs) }; }
        mp_ss operator*(const double &rhs) { return { value * rhs }; }

        mp_ss operator/(const int &rhs)    { return { value / static_cast<double>(rhs) }; }
        mp_ss operator/(const double &rhs) { return { value / rhs }; }
    
      // Constructors from different classes. Allows std::array<Tesla, n> a = { 1, 2, 3, ... }
        mp_ss()                     {this->value = 0;}
        mp_ss(const double &doub)   {this->value = doub;}
        mp_ss(const int &iint)      {this->value = static_cast<double>(iint);}
};
/*********** VELOCTI SPATIAL DERIVATIVE 1/S **********/

class _ps{
    public:

        double value;

      // Conversion operator to allow [class] a = mps b
        operator double()       {return this->value;}

      // Assignement operator to allow mps a = [class] b
        void operator=(const double &rhs)       {this->value = rhs;}
        void operator=(const int &rhs)          {this->value = static_cast<double>(rhs);}

      // MDAS
        _ps operator+(const _ps &rhs) const    { return { value + rhs.value }; }
        _ps operator-(const _ps &rhs) const    { return { value - rhs.value }; }
        _ps operator*(const _ps &rhs) const    { return { value * rhs.value }; }
        _ps operator/(const _ps &rhs) const    { return { value / rhs.value }; }

        _ps operator+(const int &rhs)    { return { value + static_cast<double>(rhs) }; }
        _ps operator+(const double &rhs) { return { value + rhs }; }

        _ps operator-(const int &rhs)    { return { value - static_cast<double>(rhs) }; }
        _ps operator-(const double &rhs) { return { value - rhs }; }

        _ps operator*(const int &rhs)    { return { value * static_cast<double>(rhs) }; }
        _ps operator*(const double &rhs) { return { value * rhs }; }

        _ps operator/(const int &rhs)    { return { value / static_cast<double>(rhs) }; }
        _ps operator/(const double &rhs) { return { value / rhs }; }
    
      // Constructors from different classes. Allows std::array<Tesla, n> a = { 1, 2, 3, ... }
        _ps()                     {this->value = 0;}
        _ps(const double &doub)   {this->value = doub;}
        _ps(const int &iint)      {this->value = static_cast<double>(iint);}
};

/********** TIME **********/

class Seconds{
    public:

        double value;

      // Conversion operator to allow [class] a = Seconds b
        operator double()       {return this->value;}
        operator Years()        {return this->value * units::stoY;}

      // Assignement operator to allow Seconds a = [class] b
        void operator=(const double &rhs)       {this->value = rhs;}
        void operator=(const int &rhs)          {this->value = static_cast<double>(rhs);}
        void operator=(const Years &rhs)        {this->value = units::Ytos * rhs.value;}

      // MDAS
        Seconds operator+(const Seconds &rhs) const    { return { value + rhs.value }; }
        Seconds operator-(const Seconds &rhs) const    { return { value - rhs.value }; }
        Seconds operator*(const Seconds &rhs) const    { return { value * rhs.value }; }
        Seconds operator/(const Seconds &rhs) const    { return { value / rhs.value }; }

        Seconds operator+(const int &rhs)    { return { value + static_cast<double>(rhs) }; }
        Seconds operator+(const double &rhs) { return { value + rhs }; }
        Seconds operator+(const Years &rhs)  { return { value + rhs.value * units::Ytos}; }

        Seconds operator-(const int &rhs)    { return { value - static_cast<double>(rhs) }; }
        Seconds operator-(const double &rhs) { return { value - rhs }; }
        Seconds operator-(const Years &rhs)  { return { value - rhs.value * units::Ytos}; }

        Seconds operator*(const int &rhs)    { return { value * static_cast<double>(rhs) }; }
        Seconds operator*(const double &rhs) { return { value * rhs }; }
        Seconds operator*(const Years &rhs)  { return { value * rhs.value * units::Ytos}; }

        Seconds operator/(const int &rhs)    { return { value / static_cast<double>(rhs) }; }
        Seconds operator/(const double &rhs) { return { value / rhs }; }
        Seconds operator/(const Years &rhs)  { return { value / (rhs.value * units::Ytos)}; }

      // Constructors from different classes. Allows std::array<Tesla, n> a = { 1, 2, 3, ... }
        Seconds()                     {this->value = 0;}
        Seconds(const double &doub)   {this->value = doub;}
        Seconds(const int &iint)      {this->value = static_cast<double>(iint);}
        Seconds(const Years &Y)       {this->value = units::Ytos * Y.value;}
        Seconds(const Seconds &s)     {this->value = s.value;}
};

class Years{
    public:

        double value;

      // Conversion operator to allow [class] a = Tesla b
        operator double()       {return this->value;}
        operator Seconds()      {return this->value * units::Ytos;}

      // Assignement operator to allow Tesla a = [class] b
        void operator=(const double &rhs)       {this->value = rhs;}
        void operator=(const int &rhs)          {this->value = static_cast<double>(rhs);}
        void operator=(const Seconds &rhs)      {this->value = units::stoY * rhs.value;}

      // MDAS
        Years operator+(const Years &rhs) const    { return { value + rhs.value }; }
        Years operator-(const Years &rhs) const    { return { value - rhs.value }; }
        Years operator*(const Years &rhs) const    { return { value * rhs.value }; }
        Years operator/(const Years &rhs) const    { return { value / rhs.value }; }

        Years operator+(const int &rhs)    { return { value + static_cast<double>(rhs) }; }
        Years operator+(const double &rhs) { return { value + rhs }; }
        Years operator+(const Seconds &rhs)  { return { value + rhs.value * units::stoY}; }

        Years operator-(const int &rhs)    { return { value - static_cast<double>(rhs) }; }
        Years operator-(const double &rhs) { return { value - rhs }; }
        Years operator-(const Seconds &rhs)  { return { value - rhs.value * units::stoY}; }

        Years operator*(const int &rhs)    { return { value * static_cast<double>(rhs) }; }
        Years operator*(const double &rhs) { return { value * rhs }; }
        Years operator*(const Seconds &rhs)  { return { value * rhs.value * units::stoY}; }

        Years operator/(const int &rhs)    { return { value / static_cast<double>(rhs) }; }
        Years operator/(const double &rhs) { return { value / rhs }; }
        Years operator/(const Seconds &rhs)  { return { value / (rhs.value * units::stoY)}; }

      // Constructors from different classes. Allows std::array<Tesla, n> a = { 1, 2, 3, ... }
        Years()                     {this->value = 0;}
        Years(const double &doub)   {this->value = doub;}
        Years(const int &iint)      {this->value = static_cast<double>(iint);}
        Years(const Seconds &s)     {this->value = units::stoY * s.value;}
        Years(const Years &Y)       {this->value = Y.value;}
};

/********** ENERGY **********/

class eV{
    public:

        double value;

      // Conversion operator to allow [class] a = mps b
        operator double()       {return this->value;}

      // Assignement operator to allow mps a = [class] b
        void operator=(const double &rhs)       {this->value = rhs;}
        void operator=(const int &rhs)          {this->value = static_cast<double>(rhs);}

      // MDAS
        eV operator+(const eV &rhs) const    { return { value + rhs.value }; }
        eV operator-(const eV &rhs) const    { return { value - rhs.value }; }
        eV operator*(const eV &rhs) const    { return { value * rhs.value }; }
        eV operator/(const eV &rhs) const    { return { value / rhs.value }; }

        eV operator+(const int &rhs)    { return { value + static_cast<double>(rhs) }; }
        eV operator+(const double &rhs) { return { value + rhs }; }

        eV operator-(const int &rhs)    { return { value - static_cast<double>(rhs) }; }
        eV operator-(const double &rhs) { return { value - rhs }; }

        eV operator*(const int &rhs)    { return { value * static_cast<double>(rhs) }; }
        eV operator*(const double &rhs) { return { value * rhs }; }

        eV operator/(const int &rhs)    { return { value / static_cast<double>(rhs) }; }
        eV operator/(const double &rhs) { return { value / rhs }; }
      
      // Extra MDAS
        mps operator/(const microGauss &rhs) const {
            return {value / (rhs.value * units::mGtoT)};
        }
    
      // Constructors from different classes. Allows std::array<Tesla, n> a = { 1, 2, 3, ... }
        eV()                     {this->value = 0;}
        eV(const double &doub)   {this->value = doub;}
        eV(const int &iint)      {this->value = static_cast<double>(iint);}
};



#endif