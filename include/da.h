#ifndef __DA_H__
#define __DA_H__

#include <include/csv.h>
#include <include/infra.h>
#include <include/lot.h>
#include <vector>

/**
 * @struct da_station_t
 * @brief store the information about a D/A station
 *
 * da_station_t is used to store the information such as how much die will be
 * processed in D/A station and the lot which is arrived or unarrived in D/A
 * station. These information is used to predict that if lot will be processed
 * in D/A station in 24 hours or not.
 *
 * @var fcst : the amount of die expect to be processed in D/A station in 24
 * hours
 * @var act : the actual amount of die which is processed in D/A station.
 * @var remain : the remain amount of die which is going to be processed.
 * @var upm : unit per minute. How many die are processed in a minute.
 * @var finished : does da_station_t finish arranging the production capacity.
 * @var arrived : store the lot which is waiting in D/A station
 * @var unarrived : store the lot which is unarrived to D/A station.
 * @var remaining : store the lot which will not be processed in D/A in 24
 * hours.
 */
typedef struct da_station_t {
    double fcst;
    double act;
    double remain;
    double upm;
    double time;
    bool finished;
    std::vector<lot_t> arrived;
    std::vector<lot_t> unarrived;
    std::vector<lot_t> remaining;
} da_station_t;

/**
 * @class da_stations_t
 * @brief used to manage all da_station_t objects
 *
 * da_stations_t object manages all da_station_t objects. The program use the
 * production capacity to predict whether the lot will pass D/A station in 24
 * hours or not. The classification lot in D/A station is by its recipe(bd_id).
 * If lot traverses to D/A station, lots will add in queue and
 */
class da_stations_t
{
private:
    std::map<std::string, da_station_t> _da_stations_container;

    std::vector<lot_t> daDistributeCapacity(
        da_station_t &da);  // for single da_station

    std::vector<lot_t> getSubLot(std::vector<lot_t> lots);

    std::vector<lot_t> _parent_lots;

public:
    da_stations_t(csv_t fcst, bool strict = false);

    /**
     * setFcst ()
     *
     * return 0 : -n
     */
    int setFcst(csv_t fcst, bool strict = false);

    bool addArrivedLotToDA(lot_t &lot);
    bool addUnarrivedLotToDA(lot_t &lot);

    std::vector<lot_t> distributeProductionCapacity();

    void removeAllLots();

    std::vector<lot_t> getParentLots();
};

inline std::vector<lot_t> da_stations_t::getParentLots()
{
    return _parent_lots;
}

#endif
