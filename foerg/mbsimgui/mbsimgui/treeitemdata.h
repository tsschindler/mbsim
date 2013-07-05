#ifndef _TREEITEMDATA__H_
#define _TREEITEMDATA__H_

#include <string>

class TreeItemData {
  public:
    virtual ~TreeItemData() {}
    virtual const std::string& getName() const = 0;
    virtual std::string getValue() const = 0;
    virtual void setName(const std::string &data) = 0;
    virtual void setValue(const std::string &data) {}
    virtual bool isRemovable() {return true;}
};

#endif