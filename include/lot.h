#ifndef __LOT_H__
#define __LOT_H__

#include <include/csv.h>
#include <include/entity.h>
#include <include/infra.h>
#include <include/job.h>
#include <algorithm>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

class lot_t
{
protected:
    std::string _route;
    std::string _lot_number;
    std::string _pin_package;
    std::string _recipe;
    std::string _prod_id;
    std::string _process_id;
    std::string _bom_id;
    std::string _part_id;
    std::string _part_no;
    std::string _urgent;
    std::string _customer;

    int _qty;
    int _oper;
    int _lot_size;
    int _amount_of_wires;
    int _amount_of_tools;

    bool _hold;
    bool _mvin;
    bool _is_sub_lot;

    double _queue_time;  // for all queue time;
    double _fcst_time;   // for DA fcst time
    double _outplan_time;

    bool _finish_traversal;

    std::vector<std::string> _log;
    std::vector<std::string> _can_run_models;
    std::vector<std::string> _can_run_locations;
    std::vector<std::string> _can_run_entities;

    std::map<std::string, double> _uphs;
    std::map<std::string, double> _model_process_times;
    std::map<std::string, double> _entity_process_times;

public:
    int tmp_oper;
    bool tmp_mvin;

    lot_t() {}

    /**
     * constructor for lot_t
     *
     * This constructor for lot_t has only one parameter whose type is
     * std::map<std::string, std::string>. The parameter is used to store the
     * relationship between key and data from dataframe. For example,
     * elements["route"] == "BGA321", elements["lot_number"] == "AASJSDKA01".
     * The mapping relationship is from WIP(work in process) sheet. The
     * constructor will initialize lots of data members. After initalizing the
     * data members, the constructor will checke if data member has incorrect
     * data. If data has incorrect data, the constructor will throw exception
     * and the error message will point out which information is not provided or
     * incorrect.
     *
     * @exception : std::invalid_argument
     * @param elements : mapping relationship between key and data for
     * datamember.
     */
    lot_t(std::map<std::string, std::string> elements) noexcept(false);

    /**
     * createSublots() - create sub lots
     *
     * The function is used to creat sub-lot, if lot isn't sub-lot. The number
     * of sublots is ⌈(qty / lot_size)⌉.
     *
     * @return a vector of lot_t object splited from the lot.
     */
    std::vector<lot_t> createSublots();

    /**
     * checkFormation() - check the basic formation of object is correct or not.
     * checkFormation is used to check the basic fomation of data member is
     * correct or not. Basic fomation of data member is from WIP(work in
     * process) sheet. If basic formation is incorrect or not provided, the lot
     * will not be involved in scheduling plan.
     *
     * @exception : std::invalid_argumaent
     */
    bool checkFormation();

    /**
     * addLog () - add log for this lot
     */
    void addLog(std::string _text);

    /**
     * addQueueTime () - add queue time
     *
     * The unit of queue time is minute. Lot initialize _queue_time as 0 in the
     * constructor. addQueueTime is used to add time to _queue_time. The
     * function is used in route_t::calculateQueueTime to add queue time station
     * by station.
     *
     * @param time : queue time in a station
     */
    void addQueueTime(double time);

    /**
     * setTraverseFinished () - set the traversing flag to be finished
     *
     * If the lot finish traversing the route, the traversing flag need to be
     * set by using this function. The lot will no longer traverse the route
     * again.
     */
    void setTraverseFinished();

    /**
     * setProcessId () - setup the process id
     *
     * The process id is mapped by product_id.
     *
     * @param a std::string type parameter
     */
    void setProcessId(std::string process_id);

    /**
     * setBomId () - setup Bom Id
     *
     * Bom id is mapped by process id
     *
     * @param a std::string type parameter
     */
    void setBomId(std::string bom_id);

    /**
     * setLotSize () - setup lot size
     *
     * The lot can be split out to several sub-lots. The number of qty of each
     * sublot is determined by lot size.
     *
     * @exception : std::invalid_argument if parameter less than 0
     * @param lotsize : a integer type of parameter
     */
    void setLotSize(int lotsize) noexcept(false);

    /**
     * setFcstTime () - setup forcast time
     *
     * forcast time is used to predict the queue time in D/A station
     *
     * @param time : a double type of variable which is in minute unit
     */
    void setFcstTime(double time);

    /**
     * setPartId () - setup part_id
     *
     * part id is used to determine what kind of tool should be used for this
     * lot
     *
     * @param partid : a std::string type of variable.
     */
    void setPartId(std::string partid);

    /**
     * setPartNo () - setup part_no
     *
     * part_no is used to determine what kind of wire should be used for this
     * lot
     *
     * @param part_no : a std::string type of variable
     */
    void setPartNo(std::string part_no);

    /**
     * setAmountOfTools () - setup amount of tools can be used for this lot
     *
     * @param amount : a integer type of variable
     */
    void setAmountOfTools(int amount);

    /**
     * setAmountOfWires () - setup amount of wires which can be used for this
     * lot
     *
     * @param amount : a integer type of variable
     */
    void setAmountOfWires(int amount);

    /**
     * setCanRunModel () - setup a single model
     */
    void setCanRunModel(std::string model);

    /**
     * setCanRunModels () - setup multiple models
     */
    void setCanRunModels(std::vector<std::string> models);


    int getAmountOfTools();

    int getAmountOfWires();

    std::vector<std::string> getCanRunModels();

    // /**
    //  * getAmountOfMachines () - get the amount of machines which can process
    //  * this lot
    //  */
    // int getAmountOfMachines();

    /**
     * oper () - current oper of this lot
     */
    int oper();

    /**
     * qty () - quantity of this lot
     */
    int qty();

    /**
     * hold () - Is this lot held or not
     *
     * @return true if lot is held
     */
    bool hold();

    /**
     * mvin () - Is this lot move into the machine or not
     *
     * @return true if lot is moving into the machine
     */
    bool mvin();

    /**
     * isSubLot () - Is this lot a sub-lot
     *
     * if length of lot number greater then 8 which means that the lot isn't a
     * sub-lot
     *
     * @return true if lot is a sub-lot
     */
    bool isSubLot();

    /**
     * isTraversalFinished () - Is the traversal finished or not?
     *

     * Lot will traverse each station in the route to sum the queue time and
     determine the arrival time to W/B station. isTraversalFinish is going to
     describe that if traversal is finished or not.
     *
     * @return true if traversal is finished.
     */
    bool isTraversalFinish();

    /**
     * route () - get the route of this lot
     *
     * @return std::string type data
     */
    std::string route();

    /**
     * log () - get the log of lot
     *
     * @return std::string type data
     */
    std::string log();

    /**
     * lotNumber () - get the lot_number
     *
     * lotNumber in WIP(work in process) is unique which represent an identifier
     * of a lot.
     *
     * @return std::string type data
     */
    std::string lotNumber();

    /**
     * recipe () - get the recipe of the lot
     *
     * recipe is bd_id.
     *
     * @return std::string type data
     */
    std::string recipe();

    /**
     * processId () - get the process_id of the lot
     *
     * process id is mapped by prod_id. Initially, process_id is empty because
     * this information is not in WIP. process_id is set by setProcessId.
     *
     * @return std::string type data
     */
    std::string processId();

    /**
     * prodId () - get the product_id of the lot
     *
     * product id is in WIP so that the return data will not be an empty string.
     *
     * @return std::string type data
     */
    std::string prodId();

    /**
     * bomId () - get the bom_id of the lot
     *
     * bom id, mapped by process_id, is used to determine what kind of wire
     * should be used on this lot.  Initially, bom_id is an empty string. Call
     * setBomId to setup the bom_id
     *
     * @return std::string type data
     */
    std::string bomId();

    /**
     * pin_package () - get the pin_package of the lot
     *
     * pin_package, an information of lot, is used to
     */
    std::string pin_package();

    /**
     * info () - get the infomation of the lot
     *
     * The infomation includes lot_number, route and recipe. Note that info is a
     * virtual function.
     *
     * @return std::string type data which include lot_number, route and recipe
     */
    virtual std::string info();

    /**
     * part_id () - get the part id of the lot
     *
     * @return std::string type of data.
     */
    std::string part_id();

    /**
     * part_no () - get the part number of the lot
     *
     * @return std::string type of data.
     */
    std::string part_no();

    /**
     * fcst () - get the forecast time of the lot
     *
     * @return forecast time in DA station, double data type, in minutes.
     */
    double fcst();

    /**
     * queueTime () - get the queue time of the lot
     *
     * if lot is finished traversing the route, queue time is also the arrival
     * time to W/B station.
     *
     * @return  time, double data type, in minutes
     */
    double queueTime();

    std::map<std::string, std::string> data();

    bool setUph(std::string name, double uph);

    bool setUph(csv_t uph);

    void setCanRunLocation(
        std::map<std::string, std::vector<std::string> > model_locations);

    std::vector<std::string> can_run_locations();

    bool isEntityCanRun(std::string model, std::string location);

    bool addCanRunEntity(entity_t *ent);

    std::vector<std::string> getCanRunEntities();

    job_t job();

    std::map<std::string, double> getEntitiyProcessTime();

    void clearCanRunLocation();

    std::map<std::string, double> uphs();
};

inline void lot_t::clearCanRunLocation()
{
    _can_run_locations.clear();
}

inline std::map<std::string, double> lot_t::uphs()
{
    return _uphs;
}

inline std::vector<std::string> lot_t::getCanRunEntities()
{
    return _can_run_entities;
}

inline std::vector<std::string> lot_t::can_run_locations()
{
    return _can_run_locations;
}

inline void lot_t::setBomId(std::string bom_id)
{
    _bom_id = bom_id;
}

inline std::string lot_t::part_no()
{
    return _part_no;
}

inline bool lot_t::isTraversalFinish()
{
    return _finish_traversal;
}

inline int lot_t::oper()
{
    return _oper;
}

// inline int lot_t::getAmountOfMachines()
// {
//     return std::min(_amount_of_tools, _amount_of_wires);
// }
inline int lot_t::getAmountOfTools()
{
    return _amount_of_tools;
}

inline int lot_t::getAmountOfWires()
{
    return _amount_of_wires;
}

inline std::string lot_t::route()
{
    return _route;
}


inline bool lot_t::hold()
{
    return _hold;
}

inline bool lot_t::mvin()
{
    return _mvin;
}

inline void lot_t::addLog(std::string _text)
{
    _log.push_back(_text);
}

inline std::string lot_t::log()
{
    std::string text;
    iter(_log, i) { text += _log[i] + " "; }
    return text;
}

inline std::string lot_t::lotNumber()
{
    return _lot_number;
}


inline void lot_t::setTraverseFinished()
{
    _finish_traversal = true;
}


inline std::string lot_t::recipe()
{
    return _recipe;
}

inline void lot_t::setProcessId(std::string process_id)
{
    _process_id = process_id;
}

inline std::string lot_t::processId()
{
    return _process_id;
}

inline std::string lot_t::prodId()
{
    return _prod_id;
}

inline std::string lot_t::bomId()
{
    return _bom_id;
}


inline void lot_t::setLotSize(int lotsize) noexcept(false)
{
    if (lotsize < 0)
        throw std::invalid_argument("lot size is less then 0");
    _lot_size = lotsize;
}

inline int lot_t::qty()
{
    return _qty;
}

inline bool lot_t::isSubLot()
{
    return _is_sub_lot;
}

inline void lot_t::setFcstTime(double time)
{
    _fcst_time += time;
    _queue_time += time;
}


inline void lot_t::addQueueTime(double time)
{
    _queue_time += time;
}


inline double lot_t::fcst()
{
    return _fcst_time;
}


inline double lot_t::queueTime()
{
    return _queue_time;
}

inline std::string lot_t::info()
{
    return "{Lot Number : " + _lot_number + ", route : " + _route +
           ", recipe : " + _recipe + "}";
}


inline void lot_t::setPartId(std::string partid)
{
    _part_id = partid;
}


inline void lot_t::setPartNo(std::string part_no)
{
    _part_no = part_no;
}


inline std::string lot_t::part_id()
{
    return _part_id;
}


inline void lot_t::setAmountOfTools(int tool)
{
    _amount_of_tools = tool;
}


inline void lot_t::setAmountOfWires(int amount)
{
    _amount_of_wires = amount;
}

inline void lot_t::setCanRunModel(std::string model)
{
    _uphs[model] = 0;
    _model_process_times[model] = 0;
    if (find(_can_run_models.begin(), _can_run_models.end(), model) !=
        _can_run_models.end())
        _can_run_models.push_back(model);
}

inline void lot_t::setCanRunModels(std::vector<std::string> models)
{
    iter(models, i)
    {
        _uphs[models[i]] = 0;
        _model_process_times[models[i]] = 0;
    }
    _can_run_models = models;
}

inline bool lot_t::setUph(std::string model, double uph)
{
    if (uph == 0) {
        _uphs.erase(model);
        _model_process_times.erase(model);
        return false;
    } else {
        _uphs.at(model) = uph;
        _model_process_times[model] = (_qty / uph) * 60;
    }
    return true;
}

inline std::vector<std::string> lot_t::getCanRunModels()
{
    return _can_run_models;
}

inline std::string lot_t::pin_package()
{
    return _pin_package;
}

typedef struct {
    std::string wire_tools_name;
    std::string wire_name;
    std::string tool_name;
    unsigned long lot_amount;
    int tool_amount;
    int wire_amount;
    int machine_amount;
    std::map<std::string, int> models_statistic;
    std::vector<entity_t *> entities;
    std::vector<std::string> entity_names;
    std::vector<lot_t *> lots;
} lot_group_t;

bool lot_group_comparision(lot_group_t g1, lot_group_t g2);


class lots_t
{
protected:
    std::vector<lot_t> lots;
    std::map<std::string, std::vector<lot_t *> > wire_lots;
    std::map<std::string, std::vector<lot_t *> > tool_lots;
    std::map<std::string, std::vector<lot_t *> > tool_wire_lots;
    std::map<std::string, int> amount_of_wires;
    std::map<std::string, int> amount_of_tools;

    std::map<std::string, int> sta_models(
        std::map<std::string, std::map<std::string, std::vector<entity_t *> > >
            ents);  // location -> amount
    std::map<std::string, int> sta_models(
        std::map<std::string, std::vector<entity_t *> > loc_ents);

    bool toolWireLotsHasLots();

public:
    void addLots(std::vector<lot_t> lots);

    std::vector<lot_group_t> round(entities_t machines);
    std::vector<std::vector<lot_group_t> > rounds(entities_t ents);

    inline std::map<std::string, int> amountOfWires()
    {
        return amount_of_wires;
    }
    inline std::map<std::string, int> amountOfTools()
    {
        return amount_of_tools;
    }
};

#endif
