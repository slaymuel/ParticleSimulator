#pragma once

#include <Eigen/Dense>

class Particle{

    private:
    //static int _count;

    public:
    Particle(){
        //this->index = _count;
        //_count++;
    }
    virtual ~Particle(){}

    Eigen::VectorXd pos;
    Eigen::VectorXd* chargePos;
    int index;

    virtual void translate() = 0;
    virtual void rotate() = 0;
};

class RPM : public Particle{

    public:
    RPM(){}
    ~RPM(){}

    void translate(){
    }

    void rotate(){
    }
};

class ARPM : public Particle{
    public:
    ARPM(){}
    ~ARPM(){}

    void translate(){
    }

    void rotate(){
    }
};