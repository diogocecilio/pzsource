#ifndef TPZBNDCONDBASE_H
#define TPZBNDCONDBASE_H
#include "TPZMaterialT.h"
#include "TPZBndCondT.h"



/**
 * @brief This class ensures that the boundary condition material is compatible
 * with any given interfaces.
 * @ingroup matsinglespace
 */
template<class TVar, class...Interfaces>
class TPZBndCondBase :
    public TPZMaterialT<TVar>,
    public TPZBndCondT<TVar>,
    public Interfaces... {
 protected:
    //! Default constructor
    TPZBndCondBase() = default;
    bool IsMatImpl() final{return true;}
 public:
    /*no need to specify destructor, copy/move assignment/destructors.
     the compiler can do it for us*/

    //! Constructor from an existing material setting `fVal1` and `fVal2`
    TPZBndCondBase(TPZMaterial * material,
                   int id,int type,
                   const TPZFMatrix<TVar> &val1,
                   const TPZVec<TVar> &val2);
    
    [[nodiscard]] int ClassId() const override;
    void Read(TPZStream &buf, void*context) override;
    void Write(TPZStream &buf, int withclassid) const override;

    void SetMaterial(TPZMaterial *) final;
    
    [[nodiscard]] int Dimension() const final
   // {return this->fMaterial->Dimension();}
    {return this->fMaterial->Dimension() - 1;}
    [[nodiscard]] int NStateVariables() const final
    {return this->fMaterial->NStateVariables();}

    [[nodiscard]] int Id() const override {
        return TPZMaterial::Id();
    }

       [[nodiscard]] int NSolutionVariables(int var) const final
    { return this->fMaterial->NSolutionVariablesBC(var);}
    
    [[nodiscard]] int VariableIndex(const std::string &name) const final
     { return this->fMaterial->VariableIndex(name);}
     
    [[nodiscard]] TPZMaterial* NewMaterial() const override {
        return new TPZBndCondBase(*this);
    }
};

template<class TVar, class...Interfaces>
TPZBndCondBase<TVar, Interfaces...>::TPZBndCondBase(TPZMaterial * material,
                                                    int id,int type,
                                                    const TPZFMatrix<TVar> &val1,
                                                    const TPZVec<TVar> &val2) :
    TPZMaterialT<TVar>(id), TPZBndCondT<TVar>(type,val1,val2)
{
    SetMaterial(material);
}

template<class TVar, class...Interfaces>
void TPZBndCondBase<TVar, Interfaces...>::SetMaterial(TPZMaterial *mat){
    auto *tmp =
        dynamic_cast<TPZMaterialT<TVar>*>(mat);
    if(!tmp){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nIncompatible material type. Aborting...\n";
        DebugStop();
    }
    this->fMaterial = tmp;
    // the next statement will call the method SetMaterialImpl(mat) for all interfaces
    // it is called a fold expression (from c++17 on)
    // https://en.cppreference.com/w/cpp/language/fold
    (Interfaces::SetMaterialImpl(mat),...);
}


template<class TVar, class...Interfaces>
int TPZBndCondBase<TVar, Interfaces...>::ClassId() const{
    constexpr int nInterfaces = sizeof...(Interfaces);
    int id = Hash("TPZBndCondBase") ^
        TPZMaterialT<TVar>::ClassId() << 1 ^
        TPZBndCondT<TVar>::ClassId() << 2;
    
//    if constexpr (nInterfaces>0){
//        TPZVec<int> classIds = {Interfaces::ClassId()...};
//        for(auto count =0;count< nInterfaces;count++){
//            id = id ^ (classIds[count] << (count+3));
//        }
//    }
    return id;
}



template<class TVar, class...Interfaces>
void TPZBndCondBase<TVar, Interfaces...>::Read(
    TPZStream& buf, void* context)
{
    TPZMaterialT<TVar>::Read(buf,context);
    TPZBndCondT<TVar>::Read(buf,context);
    /* The following will perform calls to all the Read methods of the interfaces. 
       This is a c++17 addition called fold expressions.*/
//    (Interfaces::Read(buf,context),...);
}

template<class TVar, class...Interfaces>
void TPZBndCondBase<TVar, Interfaces...>::Write(
    TPZStream& buf, int withclassid) const
{
    TPZMaterialT<TVar>::Write(buf,withclassid);
    TPZBndCondT<TVar>::Write(buf,withclassid);
    /* The following will perform calls to all the Read methods of the interfaces. 
       This is a c++17 addition called fold expressions.*/
//    (Interfaces::Write(buf,withclassid),...);
}
#endif
