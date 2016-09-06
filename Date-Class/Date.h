//
//  Date.h
//  Financial Computing in C++
//
//  Created by Huachen Qin on 10/25/15.
//  Copyright © 2015 Huachen Qin. All rights reserved.
//

#ifndef Date_h
#define Date_h

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

#endif /* Date_h */
