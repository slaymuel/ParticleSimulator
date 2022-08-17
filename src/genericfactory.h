#pragma once

namespace Simulator {

template<class ObjectType, typename Identifier, typename Creator>
class GenericFactory{
    private:
    using CallbackMap = std::map<Identifier, Creator>;
    CallbackMap callbacks;

    public:


    bool registerObject(const Identifier& id, Creator callback){
        // value_type creates a pair of type std::pair<SamplerTypes::, Callback>
        bool added = callbacks.insert(std::make_pair(id, callback)).second;
        return added;
    }
    std::unique_ptr<ObjectType> createObject(Identifier sampler_type, std::string name, std::vector<double> args){
        typename CallbackMap::const_iterator i = callbacks.find(sampler_type);
        if(i == callbacks.end()){
            Logger::Log<Logger::LogLevel::ERROR>("Generic factory cannot find the specified object!");
            return nullptr;
        }

        return (i->second)(name, args);
    }
};

} // end of namespace Simulator