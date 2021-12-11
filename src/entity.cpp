#include <include/entity.h>
#include <include/machine.h>
#include <limits.h>
#include <sched.h>
#include <algorithm>
#include <iterator>
#include <stdexcept>

using namespace std;

machine_t to_machine(entity_t ent)
{
    return machine_t{
        .base = {.machine_no = convertEntityNameToUInt(ent.entity_name),
                 .size_of_jobs = 0,
                 .avaliable_time = (unsigned int) ent.recover_time},
        .tool = ent.tool,
        .wire = ent.wire};
}

entities_t::entities_t(string time)
{
    _time = 0;
    setTime(time);
}

void entities_t::setTime(string time)
{
    if (time.length()) {
        _time = timeConverter(time);
    }
}


void entities_t::addMachine(map<string, string> elements)
{
    if (elements["recover_time"].length() == 0) {
        throw std::invalid_argument("recover time is empty");
    }
    double recover_time = timeConverter(elements["recover_time"]);


    entity_t *ent = new entity_t();
    string model, location;
    model = elements["model"];
    location = elements["location"];
    if (ent) {
        unsigned int no = convertEntityNameToUInt(elements["entity"]);
        *ent = entity_t{.recover_time = recover_time,
                        .outplan_time = recover_time,
                        .entity_name = elements["entity"],
                        .model_name = model,
                        .location = location,
                        .hold = false,
                        .tool = NULL,
                        .wire = NULL};
        ent->name.data.number[0] = no;
        ent->name.text_size = 4;
        ent->name.number_size = 1;

        _ents.push_back(ent);

        _entities[model][location].push_back(ent);
        _loc_ents[location].push_back(ent);
        if (_model_locations.count(model) == 0) {
            _model_locations[model] = vector<string>();
        }

        if (find(_model_locations[model].begin(), _model_locations[model].end(),
                 location) == _model_locations[model].end()) {
            _model_locations[model].push_back(location);
        }

    } else {
        perror("addMachine");
    }
}

void entities_t::addMachines(csv_t _machines, csv_t _location)
{
    int mrows = _machines.nrows();
    int lrows = _location.nrows();
    map<string, string> locations;  // entity->location
    for (int i = 0; i < lrows; ++i) {
        map<string, string> elements = _location.getElements(i);
        string ent = elements["entity"];
        string loc = elements["location"];
        locations[ent] = loc;
    }

    for (int i = 0; i < mrows; ++i) {
        map<string, string> elements = _machines.getElements(i);
        try {
            elements["location"] = locations.at(elements["entity"]);
            addMachine(elements);
        } catch (std::out_of_range &e) {
            // cout<<elements["entity"]<<endl;
            elements["log"] = "can't find the location of this entity";
            _faulty_machine.push_back(elements);
        } catch (std::invalid_argument &e) {
            elements["log"] = "information is loss";
            _faulty_machine.push_back(elements);
        }
    }

    iter(_ents, i)
    {
        _ents[i]->recover_time = ((_ents[i]->recover_time - _time) > 0
                                      ? (_ents[i]->recover_time - _time)
                                      : 0) /
                                 60.0;
    }
}


std::map<std::string, std::map<std::string, std::vector<entity_t *> > >
entities_t::getEntities()
{
    return _entities;
}

std::map<std::string, std::vector<entity_t *> > entities_t::getLocEntity()
{
    return _loc_ents;
}

std::vector<entity_t *> entities_t::randomlyGetEntitiesByLocations(
    std::map<std::string, int> statistic,
    int amount)
{
    vector<entity_t *> ret;

    vector<entity_t *> pool;
    for (std::map<std::string, int>::iterator it = statistic.begin();
         it != statistic.end(); it++) {
        if (it->second == 0)
            continue;
        for (unsigned int i = 0; i < _loc_ents[it->first].size(); ++i) {
            if (!_loc_ents[it->first][i]->hold) {
                pool.push_back(_loc_ents[it->first][i]);
            }
        }
    }

    random_shuffle(pool.begin(), pool.end());

    if ((unsigned) amount < pool.size()) {
        ret = vector<entity_t *>(pool.begin(), pool.begin() + amount);
    } else {
        ret = vector<entity_t *>(pool.begin(), pool.end());
    }

    iter(ret, i) { ret[i]->hold = true; }
    return ret;
}

void entities_t::reset()
{
    iter(_ents, i) { _ents[i]->hold = false; }
}

std::vector<entity_t *> entities_t::getAllEntity()
{
    return _ents;
}

std::map<std::string, std::vector<std::string> > entities_t::getModelLocation()
{
    return _model_locations;
}
unsigned int convertEntityNameToUInt(string name)
{
    union {
        char text[4];
        unsigned int number;
    } data;
    string substr = name.substr(name.length() - 4);
    strncpy(data.text, substr.c_str(), 4);
    return data.number;
}


vector<machine_t> entities_t::machines()
{
    vector<machine_t> ms;
    iter(_ents, i)
    {
        machine_t m = machine_t{
            .base = {.machine_no =
                         convertEntityNameToUInt(_ents[i]->entity_name),
                     .size_of_jobs = 0,
                     .avaliable_time = (unsigned int) _ents[i]->recover_time},
            .tool = _ents[i]->tool,
            .wire = _ents[i]->wire};
        ms.push_back(m);
    }
    return ms;
}

std::map<std::string, machine_t *> machines_t::getMachines()
{
    return _machines;
}

ancillary_resources_t::ancillary_resources_t(std::map<std::string, int> data)
{
    tool_t *t;
    for (std::map<std::string, int>::iterator it = data.begin();
         it != data.end(); it++) {
        for (int i = 0; i < it->second; ++i) {
            t = new tool_t;
            t->time = 0;
            t->machine_no = 0;
            _tools[it->first].push_back(t);
        }
    }
}

std::vector<tool_t *> ancillary_resources_t::aRound(
    std::map<std::string, int> amounts)
{
    std::vector<tool_t *> ts;
    for (std::map<std::string, int>::iterator it = amounts.begin();
         it != amounts.end(); it++) {
        std::vector<tool_t *> tmp = aRound(it->first, it->second);
        ts += tmp;
    }

    return ts;
}

std::vector<tool_t *> ancillary_resources_t::aRound(std::string name,
                                                    int amount)
{
    sort(_tools[name].begin(), _tools[name].end(), ares_ptr_comp);
    std::vector<tool_t *> ts(_tools[name].begin(),
                             _tools[name].begin() + amount);
    return ts;
}



void machines_t::addMachines(std::vector<entity_t *> ents)
{
    iter(ents, i)
    {
        machine_t m = to_machine(*ents[i]);
        machine_t *m_ptr = new machine_t;
        *m_ptr = m;
        // FIXME : the base time isn't clear
        // the base time of the machine is 0
        m_ptr->base.avaliable_time = 0;
        _machines[ents[i]->entity_name] = m_ptr;
    }
}

std::string convertUIntToEntityName(unsigned int mno)
{
    std::string text = "B";
    union {
        char text[4];
        unsigned int number;
    } data;

    data.number = mno;
    text += data.text;
    return text;
}
