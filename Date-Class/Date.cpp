//
//  Date.cpp
//  Financial Computing in C++
//
//  Created by Huachen Qin on 10/25/15.
//  Copyright © 2015 Huachen Qin. All rights reserved.
//

#include "Date.h"
#include <iostream>

using namespace std;

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

