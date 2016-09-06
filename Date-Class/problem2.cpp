// HW7 problem 2
// Huachen Qin

#include <iostream>

using namespace std;

int min(int a, int b, int c);

class Date {
    
private:
    //Static word makes them global variables for the class, created only once
    //const makes sure it is never modified
    static const int month_nb;												//number of months
    
    static const unsigned int month_days[];								//days in each month, non-leap year
    
    static const unsigned char month_daysLeap[];							//days in each month, non-leap year
    
    
protected:																	//same as private for this class as there is no inheritance (derived class)
    static bool testLeap(const int year_);									//Test if it's a leap year
    bool isLeap;															//bool is type with only two values true or false
    int year;
    int month;
    int day;
    
public:
    
    Date();																	//Default Constructor
    Date(const int month_,const int day_,const int year_);					//Overloaded Constructor
    
    int get_month() const {return month;}									//Accessig month/day/year
    int get_day() const {return day;}
    int get_year() const {return year;}
    
    void set_month(int m){month = m;}										//Setting month/day/year
    void set_day(int d){day = d;}											//make sure you check, if the date is still consistent,
    void set_year(int y){year = y;}											//after it was modified
    
    
    
    void add_months(int m);													//add month. Different conventions are possible, e.g. 1/31/12 + 1month = 2/28/12,
    //i.e. always the end of the month. You can implment any convention you prefer.
    //Make sure to comment.
    void add_days(int m);													//add days
    
    bool testLeap() const;													//Test for Leap Year
    bool testValid() const;													//Test for valid date
    int operator-(const Date & date) const;									//Operator "-" overloading
    
    
    //The following 4 functions are static, i.e. they do not they can be access not through an object of Date class.
    //Furthermore, they cannot have the this pointer, and so cannot access the non-static data members of the Date class
    //e.g. month, year, day. Hence they require two Date objects
    
    static int Julian(const Date & date1,const Date & date2);				//Difference in days
    static double CountACT_365(const Date & date1,const Date & date2);		//Difference in years, by computing difference in days,
    //and dividing by 365 days-per-year.
    static int Count30_360L(const Date & date1,const Date & date2);			//Difference in days, assuming 30 days per month
    //and 360 days per year.
    static double Count30_360(const Date & date1,const Date & date2);		//Difference in years, by dividing Count30_360L
    //by 360 days per year.
    
    //As oppposed to the following 4 regular member functions, that have access to Date class data members,
    //hence only require one extra Date object to perform the operations.
    
    int Julian(const Date & date) const;									//same as above just called differently, and are no static member functions
    double CountACT_365(const Date & date) const;
    int Count30_360L(const Date & date) const;
    double Count30_360(const Date & date) const;
    
    
    friend ostream& operator<<(ostream& os,const Date& d);					//Stream "<<" overloading
    
    //~Date() {}
};

bool operator==(const Date& d1, const Date& d2);						//Operator "==" overloading

const unsigned int Date::month_days[] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
const unsigned char Date::month_daysLeap[] = {0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

Date::Date()
{
    month = day = year = 1;
}

Date::Date(const int month_,const int day_,const int year_)
{
    month = month_;
    day = day_;
    year = year_;
    isLeap = testLeap();
}

void Date::add_months(int m)
{
    for(int i=1;i<=m;i++)                        // increase month by value 1 at each loop
    {
        month += 1;
        if(month>12) { year++; month = month % 12; }             // if added month exceeds one year, the year will increment by one
        if(month == 2 || month == 3)
        {
            if(testLeap()) day = (day + month_daysLeap[month]) % (month_daysLeap[month-1]);   // also consider leap year and non leap year
            else day = (day + month_days[month]) % (month_days[month-1]);
        }
        else
            day = (day + month_days[month]) % (month_days[month-1]);               // the day changed according to the number of the days in corresponding month
        
    }
}

void Date::add_days(int d)
{
    for(int i=d;i>0;i--)
    {
        day += 1;
        if(testLeap())
        {
            if(day>month_daysLeap[month])
            {day = day % month_daysLeap[month]; month++;}
        }
        else
        {
            if(day>month_daysLeap[month])
            {day = day % month_days[month]; month++;}
        }
        if(month>12) {month = month % 12;year++;}
    }
}

bool Date::testLeap() const
{
    if((year % 4 == 0 && year % 100 != 0) || ( year % 400 == 0))
        return true;
    else return false;
}

bool Date::testLeap(const int yr)
{
    if((yr % 4 == 0 && yr % 100 != 0) || ( yr % 400 == 0))
        return true;
    else return false;
}

bool Date::testValid() const
{
    if(testLeap())
    {
        if(min(year,month,day)>0 && month<=12 && day<=month_daysLeap[month]) return true;
        else { cout<<"Then date is not valid"; return false; }
    }
    else
    {
        if(min(year,month,day)>0 && month<=12 && day<=month_days[month]) return true;
        else { cout<<"Then date is not valid"; return false; }
    }
}

int min(int a, int b, int c)
{
    int temp = a>b?b:a;
    return temp>c?c:temp;
}

int Date::Julian(const Date &date1, const Date &date2)
{
    int temp1 = date1.day + (153 * date1.month + 2) / 5 + 365 * date1.year + date1.year / 4 - date1.year / 100 + date1.year / 400 -32045;
    int temp2 = date2.day + (153 * date2.month + 2) / 5 + 365 * date2.year + date2.year / 4 - date2.year / 100 + date2.year / 400 -32045;
    return abs(temp1-temp2);
}

double Date::CountACT_365(const Date & date1,const Date & date2)
{
    double temp = Julian(date1, date2);
    return temp / 365;
}

int Date::Count30_360L(const Date & date1,const Date & date2)
{
    int temp1 = date1.year * 360 + date1.month * 30 + date1.day;
    int temp2 = date2.year * 360 + date2.month * 30 + date2.day;
    return abs(temp1-temp2);
}

double Date::Count30_360(const Date & date1,const Date & date2)
{
    double temp = Count30_360L(date1, date2);
    return temp / 360;
}

int Date::Julian(const Date &date) const
{
    int temp1 = date.day + (153 * date.month + 2) / 5 + 365 * date.year + date.year / 4 - date.year / 100 + date.year / 400 -32045;
    int temp2 = day + (153 * month + 2) / 5 + 365 * year + year / 4 - year / 100 + year / 400 -32045;
    return abs(temp1-temp2);
}

double Date::CountACT_365(const Date & date) const
{
    double temp = Julian(date);
    return temp / 365;
}

int Date::Count30_360L(const Date & date) const
{
    int temp1 = date.year * 360 + date.month * 30 + date.day;
    int temp2 = year * 360 + month * 30 + day;
    return temp1 - temp2;
}

double Date::Count30_360(const Date & date) const
{
    double temp = Count30_360L(date);
    return temp / 360;
}

int Date::operator-(const Date & date) const
{
    if(this == &date) return 0;
    else return Julian(*this, date);
}

bool operator==(const Date& d1, const Date& d2)
{
    if(Date::Julian(d1, d2) ==0) return true;
    else return false;
}



int main()
{
    cout<< "First create two Date class object: date1,date2"<<endl;
    Date date1(10,24,2015),date2(8,14,2012);
    cout<<"date1 is: "<<date1.get_month()<<"/"<<date1.get_day()<<"/"<<date1.get_year()<<endl;
    cout<<"date2 is: "<<date2.get_month()<<"/"<<date2.get_day()<<"/"<<date2.get_year()<<endl;
    
    cout<<endl;
    
    cout<<"Then we try to add 2 months to date1, and 30 days to date2:"<<endl;
    date1.add_months(2); date2.add_days(30);
    cout<<"The result for date1 is: "<<date1.get_month()<<"/"<<date1.get_day()<<"/"<<date1.get_year()<<endl;
    cout<<"The result for date2 is: "<<date2.get_month()<<"/"<<date2.get_day()<<"/"<<date2.get_year()<<endl;
    
    cout<<endl;
    
    cout<<"Next we try to determine if date1 and date2 are leap years"<<endl;
    if(date1.testLeap() && date1.testValid()) cout<<"date1 is a leap year"<<endl;
    else cout<<"date1 is not a leap year"<<endl;
    if(date2.testLeap() && date1.testValid()) cout<<"date2 is a leap year"<<endl;
    else cout<<"date2 is not a leap year"<<endl;
    
    cout<<endl;
    
    cout<<"Next we try to calculate the difference between date1 and date2 by two methods: "<<endl;
    cout<<"By Julian method, the difference is: "<<Date::Julian(date1, date2)<<" days, or "<<Date::CountACT_365(date1, date2)<<" years"<<endl;
    cout<<"By 360 Day Count method, the difference is: "<<Date::Count30_360L(date1, date2)<<" days, or "<<Date::Count30_360(date1, date2)<<" years"<<endl;
    
    cout<<endl;
    
    cout<<"Now we try the new defined operator - and =="<<endl;
    cout<<"date1 - date2 = "<<date1 - date2<<endl;
    if(date1 == date2) cout<<"date1 == date2 True!"<<endl;
    else cout<<"date1 == date2 False!"<<endl;
    
    cout<<endl;
    
    cout<<"End of program! Please change the coefficient of date1 and date2 for more details."<<endl;
    
}
