#include <dirent.h>
#include <include/condition_card.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <iterator>
#include <stdexcept>
#include <string>

using namespace std;

condition_cards_h::condition_cards_h() {}

condition_cards_h::condition_cards_h(int n, ...)
{
    char *ptr;
    va_list list;
    va_start(list, n);
    for (int i = 0; i < n; ++i) {
        ptr = va_arg(list, char *);
        ptr = strdup(ptr);
        formated(ptr);
        // stringToLower(ptr);
        _model_set.insert(ptr);
        free(ptr);
    }
    va_end(list);
}

std::string condition_cards_h::formated(std::string _text)
{
    // to lower and replace ' ' with '-'
    for (unsigned int i = 0; i < _text.length(); ++i) {
        if (_text[i] == ' ')
            _text[i] = '-';
        else if (_text[i] & 0x40) {
            _text[i] &= 0xDF;
        }
    }
    return _text;
}

char *condition_cards_h::formated(char *text)
{
    char *ptr = text;
    for (; *ptr; ++ptr) {
        if (*ptr == ' ')
            *ptr = '-';
        else if (*ptr & 0x40) {
            *ptr &= 0xDF;
        }
    }
    return text;
}

void condition_cards_h::addMapping(std::string std_model_name, int n, ...)
{
    va_list variables;
    va_start(variables, n);
    char *ptr;
    char *model_name = strdup(std_model_name.c_str());
    stringToLower(model_name);
    for (int i = 0; i < n; ++i) {
        ptr = va_arg(variables, char *);
        ptr = strdup(ptr);
        _model_mapping[model_name].push_back(formated(ptr));
        free(ptr);
    }
    va_end(variables);
    free(model_name);
}

string condition_cards_h::dirName(string dir)
{
    return (dir.back() == '/') ? dir : dir + '/';
}

bool isCsvFile(string filename)
{
    return filename.rfind(".csv") != string::npos;
}

bool isExcelFile(string filename)
{
    return filename.rfind(".xls") != string::npos;
}

void condition_cards_h::readConditionCardsDir(string dir_name, bool replace)
{
    dir_name = dirName(dir_name);
    DIR *dir;
    struct dirent *dent;
    struct stat file;
    string sub_dir_name;

    dir = opendir(dir_name.c_str());
    while ((dent = readdir(dir)) != NULL) {
        sub_dir_name = dir_name + dent->d_name;
        stat(sub_dir_name.c_str(), &file);
        if (S_ISDIR(file.st_mode)) {
            readConditionCards(sub_dir_name, replace);
        }
    }
}


void condition_cards_h::readConditionCards(string sub_dir_name, bool replace)
{
    DIR *dir;
    struct dirent *dent;
    struct stat file;
    card_t tmep;
    vector<string> text;
    sub_dir_name = dirName(sub_dir_name);
    dir = opendir(sub_dir_name.c_str());
    string path;

    while ((dent = readdir(dir)) != NULL) {  // check if dent is file or dir
        path = sub_dir_name + dent->d_name;
        stat(path.c_str(), &file);
        if (S_ISREG(file.st_mode)) {
            if (isCsvFile(path)) {
                text = split(dent->d_name, '-');
                card_t card =
                    readConditionCard(path, text.at(text.size() - 3),
                                      stoul(text.at(text.size() - 2)));

                // check if card.recipe exist, if yes, check card.oper exiset
                if (_models.count(card.recipe) > 0 &&
                    _models[card.recipe].count(card.oper) != 0) {
                    if (replace) {
                        _models[card.recipe][card.oper] = card;
                    } else {
                        continue;
                    }
                } else {  // recipe and oper are inexist
                    _models[card.recipe][card.oper] = card;
                }
            } else if (isExcelFile(path)) {
                // do something
            }
        }
        // printf("%s\n", path.c_str());
    }
}

card_t condition_cards_h::readConditionCard(string filename,
                                            string recipe,
                                            unsigned int oper)
{
    csv_t card(filename, "r", true, false);
    vector<vector<string> > data = card.getData();
    int idx = -1;
    iter(data, i)
    {
        if (data[i][0].compare("Capillary life time") == 0) {
            idx = i;
            break;
        }
    }

    if (idx < 0) {
        string err_msg = "Cannot find models in the file " + filename;
        throw out_of_range(err_msg);
    }

    idx += 1;  // search from the next line of "Capillary life time"

    vector<string> models = card.getRow(idx);
    vector<string> result = card.getRow(idx + 1);
    vector<string> can_run_models;
    string temp;
    char *temp_model;
    unsigned int range = min(models.size(), result.size());
    range = (13 > range) ? range : 13;

    try {
        if (_model_set.empty()) {
            iter_range(result, i, 1, range)
            {
                if (models.at(i).length() && result.at(i).length()) {
                    // check if
                    can_run_models.push_back(models.at(i));
                }
            }
        } else {
            iter(result, i)
            {
                temp = models.at(i);
                temp_model = strdup(temp.c_str());
                stringToLower(temp_model);
                // check _model_mapping
                if (_model_mapping.count(temp_model) > 0 &&
                    result.at(i).length()) {
                    can_run_models += _model_mapping[temp_model];
                } else {
                    formated(temp_model);
                    if (_model_set.count(temp_model) && result.at(i).length()) {
                        can_run_models.push_back(temp_model);
                    }
                    free(temp_model);
                }
            }
        }

    } catch (std::out_of_range &e) {
    }

    return card_t{.oper = oper, .recipe = recipe, .models = can_run_models};
}

void condition_cards_h::readBdIdModelsMappingFile(string filename)
{
    csv_t bd_model(filename, "r", true, true);
    bd_model.trim(" ");
    // map<string, map<int, card_t > > models;
    map<string, string> elements;
    unsigned int nrows = bd_model.nrows();
    unsigned int oper;
    string bd_id;
    for (unsigned int i = 0; i < nrows; ++i) {
        elements = bd_model.getElements(i);
        oper = stoul(elements["oper"]);
        bd_id = elements["bd_id"];
        if (_models[bd_id].count(oper) == 0) {
            _models[bd_id][oper] = card_t{.oper = oper, .recipe = bd_id};
        }
        _models[bd_id][oper].models.push_back(formated(elements["model"]));
    }
}
