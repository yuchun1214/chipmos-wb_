#ifndef __INFRA_H__
#define __INFRA_H__

#include <cstdio>
#include <ctime>
#include <string>
#include <vector>

#include <include/def.h>  

#define iter(vec, id) for (unsigned int id = 0; id < vec.size(); ++id)


#define iter_range(vec, id, start, end) \
    for (unsigned int id = start; id < end; ++id)

std::vector<std::string> split(char *text, char delimiter);

std::string join(std::vector<std::string> strings, std::string delimiter);

template <class T>
std::vector<T> operator+(std::vector<T> &op1, std::vector<T> op2)
{
    std::vector<T> result(op1.begin(), op1.end());
    iter(op2, i) { result.push_back(op2[i]); }
    return result;
}

template <class T>
std::vector<T> operator+=(std::vector<T> &op1, std::vector<T> op2)
{
    iter(op2, i) { op1.push_back(op2[i]); }
    return op1;
}

/**
 * timeConverter() - convert text to time_t
 */
time_t timeConverter(std::string text);

void stringToLower(char *text);
void stringToUpper(char *text);

struct __info_t{
    union{
        char text[32];
        unsigned int number[8];
    }data;
    unsigned int text_size : 5;
    unsigned int number_size : 3;
};

struct __info_t to_info(std::string s);

__qualifier__ bool isSameInfo(struct __info_t info1, struct __info_t info2);


void random(double *genes, int size);

int random_range(int start, int end, int different_num);

double randomDouble();

#endif
