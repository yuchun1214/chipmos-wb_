#include <include/infra.h>
#include <cstdlib>
#include <cstring>


std::vector<std::string> split(char *text, char delimiter)
{
    char *iter = text, *prev = text;
    std::vector<std::string> data;
    while (*iter) {
        if (*iter == delimiter ||
            *iter == '\n') {  // for unix-like newline character
            *iter = '\0';
            data.push_back(prev);
            prev = ++iter;
        } else if (*iter == '\r' &&
                   *(iter + 1) ==
                       '\n') {  // for windows newline characters '\r\n'
            *iter = '\0';
            data.push_back(prev);
            iter += 2;
            prev = iter;
        } else if (*(iter + 1) == '\0') {
            data.push_back(prev);
            ++iter;
        } else {
            ++iter;
        }
    }

    return data;
}

std::string join(std::vector<std::string> strings, std::string delimiter)
{
    if (strings.size() == 0)
        return "";
    std::string s;
    iter_range(strings, i, 0, strings.size() - 1)
    {
        s += strings[i];
        s += delimiter;
    }
    s += strings[strings.size() - 1];
    return s;
}

void stringToLower(char *text)
{
    for (; *text; ++text)
        *text |= 0x20;
}

void stringToUpper(char *text)
{
    for (; *text; ++text)
        *text ^= 0x20;
}

time_t timeConverter(std::string text)
{
    // TODO : convert text to time;
    struct tm _tm;
    sscanf(text.c_str(), "%d/%d/%d %d:%d", &_tm.tm_year, &_tm.tm_mon, &_tm.tm_mday, &_tm.tm_hour, &_tm.tm_min);
    _tm.tm_sec = 0;
    _tm.tm_isdst = false;

    return mktime(&_tm);
}

bool isSameInfo(struct __info_t info1, struct __info_t info2){
    if(info1.number_size != info2.number_size){
        return false;
    }else{
        for(unsigned int i =0; i < info1.number_size; i++)
            if(info1.data.number[i] != info2.data.number[i])
                return false;
    }
    return true;
}


void random(double *genes, int size)
{
    for (int i = 0; i < size; ++i) {
        genes[i] = (double) rand() / (double) RAND_MAX;
    }
}

int random_range(int start, int end, int different_num)
{
    if (different_num < 0) {
        return start + rand() % (end - start);
    } else {
        int rnd = start + (rand() % (end - start));
        while (rnd == different_num) {
            rnd = start + (rand() % (end - start));
        }
        return rnd;
    }
}

double randomDouble(){
    return (double) rand() / (double) RAND_MAX;
}

struct __info_t to_info(std::string s){
    struct __info_t info;
    unsigned text_size = s.length() > 32 ? 32 : s.length();
    memset(info.data.number, 0, sizeof(unsigned int) * 8);
    info.text_size = text_size;
    info.number_size = 32 - __builtin_clz(text_size >> 2) + 1;

    strncpy(info.data.text, s.c_str(), info.text_size);
    return info;
}
