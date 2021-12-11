#ifndef __ENTITY_H__
#define __ENTITY_H__

#include <time.h>
#include <map>
#include <string>
#include <vector>

#include <include/csv.h>
#include <include/machine.h>

unsigned int convertEntityNameToUInt(std::string name);
std::string convertUIntToEntityName(unsigned int);

typedef struct {
    double recover_time;
    double outplan_time;
    std::string entity_name;
    machine_info_t name;
    std::string model_name;
    std::string location;
    bool hold;
    tool_t *tool;
    wire_t *wire;
} entity_t;

machine_t to_machine(entity_t ent);

class entities_t
{
private:
    std::vector<entity_t *> _ents;

    // _entities[MODEL][AREA] is a vector of entity_t object.
    std::map<std::string, std::map<std::string, std::vector<entity_t *> > >
        _entities;

    std::map<std::string, std::vector<std::string> > _model_locations;

    std::map<std::string, std::vector<entity_t *> > _loc_ents;

    std::vector<std::map<std::string, std::string> > _faulty_machine;


    time_t _time;

public:
    /**
     * entities_t () - constructor of entities_t
     *
     * The constructor will convert @b _time to time_t type
     */
    entities_t(std::string _time);

    /**
     * addMachine() - add a machine
     *
     * @b elements is a std::map container which store the relationship between
     * header and data. For example, elements[ENTITY] == "BB211",
     * elements[MODEL] = "UTC3000" and etc... This function will convert
     * elements to entity_t object.
     */
    void addMachine(std::map<std::string, std::string> elements);

    /**
     * addMachines() - add machines from dataframe
     *
     * add machines from @csv_t type dataframe.
     */
    void addMachines(csv_t machines_csv, csv_t location_csv);

    // /**
    //  * randomlyGetEntities () - randomly get the entities by model and area
    //  */
    // std::vector<entity_t> randomlyGetEntities(std::string model_name,
    //                                           std::string area,
    //                                           int amount);

    /**
     * randomlyGetEntitiesByLocations
     */
    std::vector<entity_t *> randomlyGetEntitiesByLocations(
        std::map<std::string, int> statistic,
        int amount);

    void setTime(std::string text);

    std::map<std::string, std::map<std::string, std::vector<entity_t *> > >
    getEntities();

    std::map<std::string, std::vector<entity_t *> > getLocEntity();

    std::map<std::string, std::vector<std::string> > getModelLocation();

    void reset();

    std::vector<machine_t> machines();

    std::vector<entity_t *> getAllEntity();
};


class ancillary_resources_t
{
private:
    std::map<std::string, std::vector<ares_t *> > _tools;

public:
    ancillary_resources_t(std::map<std::string, int> tools);
    std::vector<tool_t *> aRound(std::map<std::string, int> amount);
    std::vector<tool_t *> aRound(std::string, int);
};

class machines_t
{
private:
    std::map<std::string, machine_t *> _machines;

public:
    void addMachines(std::vector<entity_t *> ents);

    std::map<std::string, machine_t *> getMachines();
};


#endif
