// 2.07a    Jun 01, 05  DisplayModel: fixed bug with duplicate beta values; don't use map! //released as 2.04a

#ifndef _DATA_DESIGN_HEADER_
#define _DATA_DESIGN_HEADER_

#include <iomanip>
#include "Data.h"

enum DesignType { designPlain=0, designInteractions=1, designUndef };

class DesignParameter {
public:
    enum DesignType DesignType() const { return m_design; }

	//ctor
    DesignParameter( enum DesignType d =designPlain ) : m_design(d) {}
private:
    enum DesignType m_design;
};

inline std::ostream& operator<<( std::ostream& o, const DesignParameter& dp ) {
    o << "DesignType: "<<( 
        dp.DesignType()==designPlain ? "Plain" 
        : dp.DesignType()==designInteractions ? "Interactions" 
        : "<undefined>" );
    return o;
}


#endif //_DATA_DESIGN_HEADER_
