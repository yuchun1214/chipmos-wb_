#include "include/time_converter.h"
#include <cstring>
#include <iostream>

#define MONTH "(0[1-9]|1[0-2]|[1-9])"
#define DATE "(0[1-9]|[1-2][0-9]|3[0-1]|[1-9])"
#define HOUR "([0-9]|[0-1][0-9]|2[0-3])"
#define MINUTE "([0-9]|[0-5][0-9])"
#define DASH "-"
#define SLASH "/"
#define COLON ":"
#define SECOND MINUTE
#define SPACE " "

using namespace std;

time_converter_base_t::time_converter_base_t(std::string _pattern)
    : pattern(_pattern), pattern_str(_pattern)
{
}

bool time_converter_base_t::isThisType(std::string text)
{
    return regex_match(text, pattern);
}

void time_converter_base_t::initialized_tm(struct tm *_tm)
{
    if (_tm) {
        memset(_tm, 0, sizeof(struct tm));
        _tm->tm_isdst = false;
    }
}

time_converter_with_dash_without_second_t::
    time_converter_with_dash_without_second_t()
    : time_converter_base_t(
          R"(\d{2})" DASH MONTH DASH DATE SPACE HOUR COLON MINUTE)
{
}

time_t time_converter_with_dash_without_second_t::operator()(string text)
{
    struct tm _tm;
    initialized_tm(&_tm);
    sscanf(text.c_str(), "%d-%d-%d %d:%d", &_tm.tm_year, &_tm.tm_mon,
           &_tm.tm_mday, &_tm.tm_hour, &_tm.tm_min);
    _tm.tm_sec = 0;
    _tm.tm_mon -= 1;
    _tm.tm_isdst = false;
    _tm.tm_year += 100;
    return mktime(&_tm);
}

time_converter_with_dash_with_second_t::time_converter_with_dash_with_second_t()
    : time_converter_base_t(
          R"(\d{2})" DASH MONTH DASH DATE SPACE HOUR COLON MINUTE COLON SECOND)
{
}

time_t time_converter_with_dash_with_second_t::operator()(string text)
{
    struct tm _tm;
    initialized_tm(&_tm);
    sscanf(text.c_str(), "%d-%d-%d %d:%d:%d", &_tm.tm_year, &_tm.tm_mon,
           &_tm.tm_mday, &_tm.tm_hour, &_tm.tm_min, &_tm.tm_sec);
    _tm.tm_mon -= 1;
    _tm.tm_isdst = false;
    _tm.tm_year += 100;
    return mktime(&_tm);
}

time_converter_with_slash_without_second_t::
    time_converter_with_slash_without_second_t()
    : time_converter_base_t(
          R"(\d{4})" SLASH MONTH SLASH DATE SPACE HOUR COLON MINUTE)
{
}

time_t time_converter_with_slash_without_second_t::operator()(string text)
{
    struct tm _tm;
    initialized_tm(&_tm);
    sscanf(text.c_str(), "%d/%d/%d %d:%d", &_tm.tm_year, &_tm.tm_mon,
           &_tm.tm_mday, &_tm.tm_hour, &_tm.tm_min);
    _tm.tm_sec = 0;
    _tm.tm_mon -= 1;
    _tm.tm_isdst = false;
    _tm.tm_year -= 1900;
    return mktime(&_tm);
}

time_converter_with_slash_with_second_t::
    time_converter_with_slash_with_second_t()
    : time_converter_base_t(R"(\d{4})" SLASH MONTH SLASH DATE SPACE HOUR COLON
                                MINUTE COLON SECOND)
{
}

time_t time_converter_with_slash_with_second_t::operator()(string text)
{
    struct tm _tm;
    initialized_tm(&_tm);
    sscanf(text.c_str(), "%d/%d/%d %d:%d:%d", &_tm.tm_year, &_tm.tm_mon,
           &_tm.tm_mday, &_tm.tm_hour, &_tm.tm_min, &_tm.tm_sec);
    _tm.tm_mon -= 1;
    _tm.tm_isdst = false;
    _tm.tm_year -= 1900;
    return mktime(&_tm);
}

time_converter_only_date_with_slash_t::time_converter_only_date_with_slash_t()
    : time_converter_base_t(R"(\d{4})" SLASH MONTH SLASH DATE)
{
}

time_t time_converter_only_date_with_slash_t::operator()(string text)
{
    struct tm _tm;
    initialized_tm(&_tm);
    sscanf(text.c_str(), "%d/%d/%d", &_tm.tm_year, &_tm.tm_mon, &_tm.tm_mday);
    _tm.tm_hour = 0;
    _tm.tm_min = 0;
    _tm.tm_sec = 0;
    _tm.tm_mon -= 1;
    _tm.tm_isdst = false;
    _tm.tm_year -= 1900;
    return mktime(&_tm);
}

time_converter_only_date_with_dash_t::time_converter_only_date_with_dash_t()
    : time_converter_base_t(R"(\d{2})" DASH MONTH DASH DATE)
{
}

time_t time_converter_only_date_with_dash_t::operator()(string text)
{
    struct tm _tm;
    initialized_tm(&_tm);
    sscanf(text.c_str(), "%d-%d-%d", &_tm.tm_year, &_tm.tm_mon, &_tm.tm_mday);
    _tm.tm_hour = 0;
    _tm.tm_min = 0;
    _tm.tm_sec = 0;
    _tm.tm_mon -= 1;
    _tm.tm_isdst = false;
    _tm.tm_year += 100;
    return mktime(&_tm);
}


vector<time_converter_base_t *> timeConverter::converters = {
    new time_converter_with_dash_without_second_t(),
    new time_converter_with_dash_with_second_t(),
    new time_converter_with_slash_with_second_t(),
    new time_converter_with_slash_without_second_t(),
    new time_converter_only_date_with_slash_t(),
    new time_converter_only_date_with_dash_t()};

time_t timeConverter::operator()(std::string text)
{
    time_t time = 0;
    for (unsigned int i = 0; i < converters.size(); ++i) {
        if (converters[i]->isThisType(text))
            time = converters[i]->operator()(text);
    }
    return time - _base_time;
}

time_t timeConverter::operator-(timeConverter tc)
{
    return tc._base_time - _base_time;
}
