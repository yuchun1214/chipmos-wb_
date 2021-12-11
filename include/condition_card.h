#ifndef __CONDITION_CARD_H__
#define __CONDITION_CARD_H__

#include <include/csv.h>
#include <include/infra.h>
#include <stdarg.h>
#include <map>
#include <set>
#include <string>
#include <type_traits>
#include <vector>

/**
 * @struct card_t
 * @brief card_t used to record the information of condtion card
 *
 * @var oper : operation
 * @var recipe : recipe for the lot
 * @var models : the machine models in the condition card
 */
typedef struct {
    unsigned int oper;
    std::string recipe;
    std::vector<std::string> models;
} card_t;

/**
 * isCsvFile () - check if filename has suffix .csv
 * @return true if filename has suffix .csv
 */
bool isCsvFile(std::string filename);

/**
 * isExcelFile () - check if filename has suffix .xls or .xlsx
 * @return true if filename has suffix .xls or .xlsx
 */
bool isExcelFile(std::string filename);

/**
 * @class condition_cards_h
 * condition_cards_h class is used to handle the information of condition cards.
 * condition_cards_h
 *
 * Condition cards are use to determine the machine models and their related
 * recipe. There are lots of condition cards in the factory so that
 * condition_cards_h has function to read a dirctory which contains bunch of
 * condition cards. The format of a model in the condition card is different.
 * Programmer can give constructor a bunch of models' names to specify the
 * text which should be recognized in the condition card so that the class
 * won't store wrong model. The name of the model is changed to lowercase and
 * space are replaced with hyphen.
 *
 */
class condition_cards_h
{
private:
    /// _models store the mapping relationship between recipe->oper->card_t
    std::map<std::string, std::map<int, card_t> > _models;

    /// _model_set stores the set of name of models which are specified in
    /// constructor for identifing the name of model in the condition cards.
    std::set<std::string> _model_set;

    /// _model_mapping store one-to-more relationship for model's name mapping
    /// For example, Maxum (Ultra) maps to Maxum and Maxum Ultra
    std::map<std::string, std::vector<std::string> > _model_mapping;

    /// store the log when executing the function member of condition_cards_h
    std::vector<std::string> _log;

    /**
     * dirName () - change the input text to directory name
     *
     * If text has / at the end, the function return text.
     * else, add / at the endn and return.
     *
     * @return a directory name
     */
    std::string dirName(std::string text);


protected:
    /**
     * addLog () - add log
     */
    virtual inline void addLog(std::string text) { _log.push_back(text); }

public:
    condition_cards_h();

    /**
     * The constructor is used to spcifiy the name of model which should be
     * identified in the condition card. The variable arguments are the name of
     * models which would be formated to be in lower case and without space.
     *
     * @param n : number of variable arguments
     * @param ... : name of models which should be identified in the condition
     * card.
     */
    condition_cards_h(int n, ...);

    /**
     * addMapping () - add an one-to-more relationship for mapping model names
     *
     * one-to-more relationship is mapping a model name to one or more model
     * names. For example, the meaning of Maxum (Ultra) is Maxum or Maxum Ultra.
     * addMapping function is used in this situation. The @b key in mapping
     * relationship passed to this function will be changed to be in lower case,
     * but not replaced space with hyphen. The data in the mapping relationship
     * will be formated to have the standard format which is lowercase and
     * without space.
     *
     * @param model_name : the key in mapping relationship
     * @param n : number of variable arguments
     * @param ... : data in the mapping relationship
     */
    void addMapping(std::string model_name, int n, ...);

    /**
     * readConditionCardsDir () - read the condition cards' directory
     *
     * The directory structure of condition card is below.
     * .
     * ├── CARD_OFFICAL
     * │   ├── AAB008YB2009A
     * │   │   └── AAB008YB2009A-2200-Criteria.csv
     * │   ├── AAB008YD2012B
     * │   │   ├── AAB008YD2012B-2200-Criteria.csv
     * │   │   └── AAB008YD2012B-3200-Criteria.csv
     * ...
     * ├── CARD_TEMP
     * │   ├── AAB024XJ6011A
     * │   │   ├── 20210413-B-AAB024XJ6011A-2200-Criteria.csv
     * │   │   ├── 20210413-B-AAB024XJ6011A-3200-Criteria.csv
     * │   │   ├── 20210413-B-AAB024XJ6011A-3400-Criteria.csv
     * │   │   └── 20210413-B-AAB024XJ6011A-3600-Criteria.csv
     * ...
     *
     * readConditionCardsDir deals with CARD_OFFICAL and CARD_TEMP. The function
     * will iterate all sub-directory and read the condition cards inside.
     * readConditionCards is in charge of reading sub-directories such as
     * AAB008Y2009, AAB008YD2012B and etc...
     *
     * @param dir_name : sub-directory name
     * @bool replace : replace the old information or not.
     */
    void readConditionCardsDir(std::string dir_name, bool replace = false);

    /**
     * readConditionCards () - read the directory which contains lots of
     * condition cards
     *
     * readConditionCards will iterate all condition cards' file name and pass
     * the path to readConditionCard function, which is in charge of reading a
     * single condition card file. The function will store data into _models
     *
     * @param sub_dir_name : name of sub-directory name
     * @param relace : replace old data if true is passed into the function
     */
    void readConditionCards(std::string sub_dir_name, bool replace = false);

    void readBdIdModelsMappingFile(std::string filename);


    /**
     * readConditionCard () - read a single condition card file
     *
     * csv_t class is in charge of I/O of condition card. The function will read
     * two rows after the row whose first elements is 'Capillary life time'
     * because there are lots of condition cards are not in good csv format.
     * Maybe there has newline character in the element string or has ^M at the
     * end of row. These situations don't influence the appearence of the file
     * in common spreadsheet software but influence in text-base file I/O. The
     * function will not use fixed row number to extract models from condition
     * card but use the 'Capillary life time' as an identifier to identify the
     * row of models.
     *
     * @param filename : the path of condition card.
     * @param recipe : the recipe of this condition card.
     * @param oper : the operation of the condition card.
     *
     * @return : card_t type object which contains the information about the
     * card, such as recipe, oper and models.
     */
    card_t readConditionCard(std::string filename,
                             std::string recipe = "",
                             unsigned int oper = 0);

    /**
     * log() - get the log from the object
     */
    std::vector<std::string> log();

    /**
     * getModels() - get the models by recipe and oper
     */
    inline card_t getModels(std::string recipe, int oper) noexcept(false)
    {
        return _models.at(recipe).at(oper);
    }

    inline std::map<std::string, std::map<int, card_t> > getModels()
    {
        return _models;
    }

    /**
     * formated () - change the format of text
     *
     * formated are used to change the format of text to be in lowercase and
     * replace the space with hyphen of the text.
     *
     * @return the formatted text
     */
    virtual std::string formated(std::string text);
    virtual char *formated(char *);
};

#endif
