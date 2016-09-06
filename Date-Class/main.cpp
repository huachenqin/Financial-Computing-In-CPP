
#include<iostream>
#include "Date.h"
using namespace std;

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
