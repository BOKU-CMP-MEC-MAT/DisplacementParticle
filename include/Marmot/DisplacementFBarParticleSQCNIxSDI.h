/* ---------------------------------------------------------------------
 * _
 * _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck,
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
 *
 * Matthias Neuner matthias.neuner@uibk.ac.at
 * Thomas Mader thomas.mader@boku.ac.at
 *
 * This file is part of the MAteRialMOdellingToolbox (marmot).
 * ---------------------------------------------------------------------
 */
#pragma once

#include "Marmot/DisplacementParticle.h"
#include "Marmot/DisplacementFBarMaterialPoint.h" // NEUER F-Bar Materialpunkt!
#include "Marmot/GenericSDIParticle.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Fastor/Fastor.h>
#include <stdexcept>

namespace Marmot::Meshfree {

  template < int nDim, int nVertices >
  class DisplacementFBarParticleSQCNIxSDI : public GenericSDIParticle< nDim, nVertices > {

  protected:
    // Shadowing: Wir definieren den MaterialPointType explizit für dieses Partikel
    using MaterialPointType = std::conditional_t< nDim == 2,
                                                  Marmot::MaterialPoints::DisplacementFBarMaterialPoint2D,
                                                  Marmot::MaterialPoints::DisplacementFBarMaterialPoint3D >;

    using TensorD  = typename GenericSDIParticle< nDim, nVertices >::TensorD;
    using TensorDD = typename GenericSDIParticle< nDim, nVertices >::TensorDD;
    
    // Vektor verwendet nun automatisch den 2D oder 3D Typ OHNE Template-Parameter <nDim>
    std::vector< std::unique_ptr< MaterialPointType > > _subdomainMaterialPoints;

    constexpr int static nStateVarsParticle = nDim * nVertices + nDim; // vertex displacements + center displacement

    double _newmark_beta;
    double _newmark_gamma;

    inline static const std::vector< std::string > _validProperties = {
      "newmark-beta beta",
      "newmark-beta gamma",
    };

  public:
    enum BodyLoadTypes { BodyForce };
    enum DistributedLoadTypes { Pressure };

    const std::unordered_map< std::string, int >& getSupportedBodyLoadTypes() const override
    {
      static const std::unordered_map< std::string, int > _supportedBodyLoadTypes = { { "BODYFORCE", BodyForce } };
      return _supportedBodyLoadTypes;
    };

    const std::unordered_map< std::string, int >& getSupportedDistributedLoadTypes() const override
    {
      static const std::unordered_map< std::string, int > _supportedDistributedLoadTypes = { { "PRESSURE", Pressure } };
      return _supportedDistributedLoadTypes;
    };

    static constexpr int nDofPerNodeU = nDim;

    virtual void setPropertyOnSubdomains( const std::string& propertyName, const double* property ) override
    {
      if ( propertyName == "newmark-beta beta" ) _newmark_beta = property[0];
      else if ( propertyName == "newmark-beta gamma" ) _newmark_gamma = property[0];
      else throw std::runtime_error( "Property " + propertyName + " not supported!" );
    };

    virtual std::vector< std::string > getSubdomainPropertyNames() const override { return _validProperties; };

    virtual int getNBaseDof() const override { return nDofPerNodeU; };

    void initializeYourselfOnSubdomains() override
    {
      for ( auto& mp : _subdomainMaterialPoints ) mp->initializeYourself();
    };

    virtual void computeBodyLoad( int type, const double* load, double* fExt, double* dExt_dQ, double timeNew, double dT ) const override
    {
      throw std::runtime_error( "Not implemented yet!" );
    }

    virtual const std::vector< std::string >& getFields() const override
    {
      static const std::vector< std::string > nodeFields = { "displacement" };
      return nodeFields;
    };

    DisplacementFBarParticleSQCNIxSDI(
      int                                                                    elementID,
      const double* nodeCoordinates,
      int                                                                    nNodeCoordiantes,
      double                                                                 volume,
      const std::string&                                                     materialName,
      const double* materialProperties,
      int                                                                    sizeMaterialProperties,
      const MarmotMeshfreeApproximation&                                     approximation,
      const typename GenericSDIParticle< nDim, nVertices >::SmoothingDomainUpdateType smoothingVolumeUpdateType );

    virtual double getVolumeUndeformed() const
    {
      double V0 = 0.0;
      for ( const auto& mp : _subdomainMaterialPoints ) V0 += mp->getVolumeUndeformed();
      return V0;
    }

    virtual double getVolumeDeformed() const
    {
      double volDeformed = 0.0;
      // Der neue MP speichert das wahre geometrische F in dY_dX() - perfekt für Volumen!
      for ( const auto& mp : _subdomainMaterialPoints ) volDeformed += mp->getVolumeUndeformed() * Fastor::determinant( mp->dY_dX() );
      return volDeformed;
    }

    virtual double getSubdomainVolume( const ParticleDomain< nDim, nVertices >& subdomain ) const override
    {
      int subdomainIndex = -1;
      for ( size_t i = 0; i < this->_subDomains.size(); i++ ) {
        if ( &subdomain == &( this->_subDomains[i] ) ) {
          subdomainIndex = static_cast< int >( i );
          break;
        }
      }
      assert( subdomainIndex != -1 && "Subdomain not found!" );
      const auto& mp = _subdomainMaterialPoints[subdomainIndex];
      return mp->getVolumeUndeformed() * Fastor::determinant( mp->dY_dX() );
    };

    virtual void acceptStateAndPositionOnSubdomains() override
    {
      for ( auto& mp : _subdomainMaterialPoints ) mp->acceptStateAndPosition();
    };

    virtual int getNumberOfRequiredStateVarsOnSubdomains() const override
    {
      int nStateVars = 0;
      for ( const auto& mp : _subdomainMaterialPoints ) nStateVars += mp->getNumberOfRequiredStateVars();
      return nStateVars;
    };

    virtual void assignStateVarsOnSubdomains( double* stateVars, int nStateVars ) override
    {
      int offset = 0;
      for ( auto& mp : _subdomainMaterialPoints ) {
        int nStateVarsSubParticle = mp->getNumberOfRequiredStateVars();
        mp->assignStateVars( stateVars + offset, nStateVarsSubParticle );
        offset += nStateVarsSubParticle;
      }
      if ( offset != nStateVars ) throw std::runtime_error( "Error: Number of state variables does not match!" );
    }

    virtual StateView getStateViewOnSubdomains( const std::string& stateName, int subdomainIndex ) const override
    {
      return _subdomainMaterialPoints[subdomainIndex]->getStateView( stateName );
    }

    virtual void computePhysicsKernelsOnSubdomains( const double* dQ, double* fInt, double* dFInt_ddQ, double timeNew, double dT ) override;
    virtual void computeDistributedLoad( int type, int surfaceID, const double* load, double* fExt, double* dExt_dQ, double timeNew, double dT ) const override;

    virtual void vci_compute_Test_P_BoundaryIntegral( double* R_AiC_RowMajor, const double* boundarySurfaceVector, int boundaryFaceID ) override
    {
      const auto [N_dAY, Y_N] = this->getIntermediateConfigurationBoundaryVector( boundaryFaceID, this->_particleDomainMain );
      Eigen::MatrixXd TBoundary = Eigen::MatrixXd::Zero( 1, this->_nNodes );

      this->_meshfreeApproximation.computeShapeFunctions( Y_N.data(), this->_assignedKernelFunctions, TBoundary.data() );

      Eigen::VectorXd PBoundary( this->_nVCIConstraints );
      Eigen::Matrix< double, nDim, 1 > Y_N_coords( Y_N.data() );
      Math::computeMonomialBasis( this->_vciOrder, Y_N_coords, PBoundary );

      for ( int A = 0; A < this->_nNodes; A++ ) {
        for ( int i = 0; i < nDim; i++ ) {
          for ( int C = 0; C < this->_nVCIConstraints; C++ ) {
            R_AiC_RowMajor[A * ( nDim * this->_nVCIConstraints ) + i * this->_nVCIConstraints + C] += 
                TBoundary( 0, A ) * PBoundary( C ) * N_dAY[i];
          }
        }
      }
    };

    virtual void getEvaluationCoordinates( double* coordinates ) const override
    {
      const int nEvalPoints = this->getNumberOfEvaluationPoints();
      Eigen::Map< Eigen::Matrix< double, nDim, Eigen::Dynamic > > coordinatesMap( coordinates, nDim, nEvalPoints );

      int i = 0;
      for ( int f = 0; f < this->_particleDomainMain.getNumberOfFaces(); f++ ) {
        const auto subcellsAttachedToFace = this->_particleDomainMain.getSubCellIndicesOnParentFace( f + 1 );
        for ( size_t sc = 0; sc < subcellsAttachedToFace.size(); sc++ ) {
          const int subcellIndex  = subcellsAttachedToFace[sc];
          coordinatesMap.col( i ) = this->_subDomains[subcellIndex].getSmoothingDomainFaceCenterCoordinates( f + 1 );
          i++;
        }
      }
    }

    virtual int getNumberOfEvaluationPoints() const override
    {
      int nEvalPoints = 0;
      for ( int i = 0; i < this->_particleDomainMain.getNumberOfFaces(); i++ ) {
        const auto subcellsAttachedToFace = this->_particleDomainMain.getSubCellIndicesOnParentFace( i + 1 );
        nEvalPoints += static_cast< int >( subcellsAttachedToFace.size() );
      }
      return nEvalPoints;
    };

    virtual void setInitialCondition( const std::string& conditionName, const double* value ) override
    {
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid initial condition" );
    };
  };

  template < int nDim, int nVertices >
  DisplacementFBarParticleSQCNIxSDI< nDim, nVertices >::DisplacementFBarParticleSQCNIxSDI(
    int                                                                    elementID,
    const double* vertexCoordinates,
    int                                                                    nVertexCoordinates,
    double                                                                 volume,
    const std::string&                                                     materialName,
    const double* materialProperties,
    int                                                                    sizeMaterialProperties,
    const Marmot::Meshfree::MarmotMeshfreeApproximation&                   approximation,
    const typename GenericSDIParticle< nDim, nVertices >::SmoothingDomainUpdateType smoothingVolumeUpdateType )
    : GenericSDIParticle< nDim, nVertices >( elementID,
                                             vertexCoordinates,
                                             nVertexCoordinates,
                                             volume,
                                             approximation,
                                             smoothingVolumeUpdateType ),
      _newmark_beta( 0. ),
      _newmark_gamma( 0. )
  {
    int materialCode = MarmotLibrary::MarmotMaterialFactory::getMaterialCodeFromName( materialName );
    MarmotMaterialSection section( materialCode, materialProperties, sizeMaterialProperties );

    for ( size_t i = 0; i < this->_subDomains.size(); i++ ) {
      const auto& initialParticleDomain = this->_subDomains[i];
      const double subV0 = initialParticleDomain.getVolumeUndeformed();
      const auto   X0    = initialParticleDomain.getCenterCoordinates();

      // unique_ptr erwartet jetzt direkt den MaterialPointType (2D oder 3D)
      _subdomainMaterialPoints.push_back( std::make_unique< MaterialPointType >( elementID, X0.data(), 1, subV0 ) );
    }

    for ( auto& mp : _subdomainMaterialPoints ) mp->assignMaterial( section );
  }

  template < int nDim, int nVertices >
  void DisplacementFBarParticleSQCNIxSDI< nDim, nVertices >::computeDistributedLoad( int           type,
                                                                                 int           boundaryFaceID,
                                                                                 const double* load_,
                                                                                 double* fExt,
                                                                                 double* dFExt_ddQ,
                                                                                 double        timeNew,
                                                                                 double        dT ) const
  {
    switch ( type ) {
    case DisplacementParticle< nDim >::Pressure: {
      constexpr int nodeBlockSize = nDim;
      const auto [N_dAY, Y_N] = this->getIntermediateConfigurationBoundaryVector( boundaryFaceID, this->_particleDomainMain );
      const auto T_Boundary = this->evaluateShapeFunctionsOnFace( this->_particleDomainMain, boundaryFaceID );
      const auto [_, dN_dY] = this->evaluateShapeFunctionsAndDerivativesForParticleDomain( this->_particleDomainMain );

      TensorD fY = N_dAY * load_[0];

      Eigen::Map< Eigen::VectorXd > P( fExt, this->_nNodes * nodeBlockSize );
      Eigen::Map< Eigen::MatrixXd > K( dFExt_ddQ, this->_nNodes * nodeBlockSize, this->_nNodes * nodeBlockSize );

      using namespace Fastor;
      using namespace FastorIndices;

      Tensor< double, nDim, nDim > Eye; Eye.eye();

      const TensorDD deltaF    = transpose( TensorDD( this->_centralDeformationGradientDelta.data() ) );
      const TensorDD deltaFInv = inverse( deltaF );
      const double   deltaJ    = determinant( deltaF );

      const TensorD f = deltaJ * transpose( deltaFInv ) % fY;
      const Tensor< double, nDim, nDim, nDim, nDim > dFInv_dF = -einsum< Ik, Ki, to_IikK >( deltaFInv, deltaFInv );
      const Tensor< double, nDim, nDim, nDim > df_dDeltaF = outer( f, transpose( deltaFInv ) ) +
                                                            deltaJ * einsum< IikK, Index< I_ > >( dFInv_dF, fY );

      TensorD r_U( 0.0 );
      for ( int A = 0; A < this->_nNodes; A++ ) {
        const int idxA_u = nodeBlockSize * A;
        r_U = T_Boundary( A ) * f;

        {
          using namespace Eigen;
          P.template segment< nDim >( idxA_u ) -= Map< Matrix< double, nDim, 1 > >( r_U.data() );
        }

        for ( int B = 0; B < this->_nNodes; B++ ) {
          const int  idxB_u  = nodeBlockSize * B;
          const auto dN_B_dY = Tensor< double, nDim >( dN_dY.col( B ).data() );
          const Tensor< double, nDim, nDim > df_ddQU_B = T_Boundary( A ) * einsum< ijk, k >( df_dDeltaF, dN_B_dY );

          {
            using namespace Eigen;
            K.template block< nDim, nDim >( idxA_u, idxB_u ) -= Map< Matrix< double, nDim, nDim > >( torowmajor( df_ddQU_B ).data() );
          }
        }
      }
      break;
    }
    default: {
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid DistributedLoad type specified" );
    }
    }
  }

  template < int nDim, int nVertices >
  void DisplacementFBarParticleSQCNIxSDI< nDim, nVertices >::computePhysicsKernelsOnSubdomains( const double* dQ,
                                                                                                double* fInt,
                                                                                                double* dFInt_ddQ,
                                                                                                double        timeNew,
                                                                                                double        dT )
  {
    using namespace Marmot::FastorIndices;
    using namespace Fastor;
    using to_jk = Fastor::OIndex< j_, k_ >;
    using ikm = Fastor::Index< i_, k_, m_ >;

    const static Tensor< double, nDim, nDim > I(
      ( Eigen::Matrix< double, nDim, nDim >() << Eigen::Matrix< double, nDim, nDim >::Identity() ).finished().data() );

    constexpr int nodeBlockSize = nDim;

    // 1. ZENTRALE KINEMATIK
    const auto [N_cen, dN_dY_cen_eigen] = this->evaluateShapeFunctionsAndDerivativesForParticleDomain( this->_particleDomainMain );

    TensorDD deltaF_cen = transpose( TensorDD(this->_centralDeformationGradientDelta.data()) );
    TensorDD inv_deltaF_cen = inverse(deltaF_cen);
    double J_cen_delta = determinant(deltaF_cen);

    for ( size_t mpNumber = 0; mpNumber < this->_subDomainShapeFunctions.size(); mpNumber++ ) {
      auto& mp = this->_subdomainMaterialPoints[mpNumber];
      auto& sd = this->_subDomainShapeFunctions[mpNumber];

      Tensor< double, nDim > du( 0.0 );
      Tensor< double, nDim, nDim > du_dY( 0.0 );

      for ( int B = 0; B < this->_nNodes; B++ ) {
        const int idxB_u = nodeBlockSize * B;
        const double N_B = sd.N( B );
        const auto dN_B_dY = Tensor< double, nDim >( sd.dN_dY.col( B ).data() );
        const auto dQU = Tensor< double, nDim >( dQ + idxB_u );

        du += N_B * dQU;
        du_dY += einsum< i, j >( dQU, dN_B_dY );
      }

      // 2. F-BAR KINEMATIK (Relative Scale)
      TensorDD deltaF_sub = I + du_dY;
      double J_sub_delta = determinant(deltaF_sub);

      // Dieser inkrementelle Skalar wird vom MP für die perfekte Total-F-Bar Integration genutzt!
      double fbarScale = std::pow( J_cen_delta / J_sub_delta, 1.0 / static_cast<double>(nDim) );

      mp->prepareYourself( timeNew, dT );
      
      // Übergebe das echte geometrische Inkrement und den F-Bar Skalar
      mp->incrementDeformation( du, du_dY, fbarScale );

      // Sende das Zentrum für den State Tracker (Optional aber sauber)
      Eigen::Matrix<double, nDim, 1> centerCoords = Eigen::Map<const Eigen::Matrix<double, nDim, 1>>(this->_particleDomainMain.getCenterCoordinates().data());
      mp->updateCenter(centerCoords, deltaF_cen);

      mp->computeYourself( timeNew, dT );

      const double density0 = mp->getDensityUndeformed();
      auto v = mp->getVelocity();
      auto a = mp->getAcceleration();

      Tensor< double, nDim, nDim > da_ddu( 0.0 );
      Marmot::TimeIntegration::newmarkBetaIntegration< nDim >( du.data(), v.data(), a.data(), dT, this->_newmark_beta, this->_newmark_gamma, da_ddu.data() );
      mp->setVelocity( v );
      mp->setAcceleration( a );

      const auto& S = mp->response.S;
      const double V0 = mp->getVolumeUndeformed();
      const auto& t = mp->tangents;

      Eigen::Map< Eigen::VectorXd > P( fInt, this->_nNodes * nodeBlockSize );
      Eigen::Map< Eigen::MatrixXd > K( dFInt_ddQ, this->_nNodes * nodeBlockSize, this->_nNodes * nodeBlockSize );

      double tr_S = 0.0;
      for(int i=0; i<nDim; ++i) tr_S += S(i,i);
      double p_tau = tr_S / static_cast<double>(nDim);

      TensorDD inv_deltaF_sub = inverse(deltaF_sub);
      
      // Hole das skalierte F-Bar Inkrement aus dem Materialpunkt!
      TensorDD deltaF_bar = mp->dx_dY_bar(); 

      // t.dS_dDeltaF ist vom MP bereits intern als dS/dDeltaF_bar skaliert worden!
      TensorDD C_F = evaluate( einsum<ijkl, kl>(t.dS_dDeltaF, deltaF_bar) );

      for ( int A = 0; A < this->_nNodes; A++ ) {

        const double T_A = sd.T( A );
        const auto dT_A_dY = Tensor< double, nDim >( sd.dT_dY.col( A ).data() );
        const Tensor< double, nDim > dT_A_dx = einsum< ji, j >( inv_deltaF_sub, dT_A_dY );

        TensorD dT_A_dY_cen = Tensor< double, nDim >( dN_dY_cen_eigen.col( A ).data() );
        TensorD dT_A_dx_cen = einsum< ji, j >( inv_deltaF_cen, dT_A_dY_cen );

        const int idxA_u = nodeBlockSize * A;

        // 3. F-BAR INTERNE KRÄFTE
        TensorD r_U_bar = einsum< i, ij >( dT_A_dx, S ) - p_tau * dT_A_dx + p_tau * dT_A_dx_cen;
        TensorD r_U = r_U_bar * V0;

        r_U += density0 * a * T_A * V0;

        {
          using namespace Eigen;
          P.template segment< nDim >( idxA_u ) += Map< Matrix< double, nDim, 1 > >( r_U.data() );
        }

        for ( int B = 0; B < this->_nNodes; B++ ) {

          const int idxB_u = nodeBlockSize * B;

          const double N_B = sd.N( B );
          const auto dN_B_dY = Tensor< double, nDim >( sd.dN_dY.col(B).data() );
          const auto dN_B_dx = evaluate( einsum< ji, j >( inv_deltaF_sub, dN_B_dY ) );

          TensorD dN_B_dY_cen = Tensor< double, nDim >( dN_dY_cen_eigen.col( B ).data() );
          TensorD dN_B_dx_cen = einsum< ji, j >( inv_deltaF_cen, dN_B_dY_cen );

          // 4. F-BAR TANGENTENMATRIX
          // Die alte Skalierung "scale_total" bzw. "scale" ist hier fbarScale
          auto dS_dqU_B = evaluate ( fbarScale * einsum < ijkl, l > ( t.dS_dDeltaF, dN_B_dY ) );

          Tensor<double, nDim, nDim, nDim> C_F_x_dN_dx = outer(C_F, dN_B_dx);
          Tensor<double, nDim, nDim, nDim> C_F_x_dN_dx_cen = outer(C_F, dN_B_dx_cen);

          auto dS_dqU_B_bar = evaluate( dS_dqU_B + (C_F_x_dN_dx_cen - C_F_x_dN_dx) / static_cast<double>(nDim) );

          Tensor<double, nDim> tr_dS_bar(0.0);
          for(int m=0; m<nDim; ++m) {
              for(int j=0; j<nDim; ++j) {
                  tr_dS_bar(m) += dS_dqU_B_bar(j, j, m);
              }
          }

          TensorDD k_UU_mat = einsum<k, ikm>( dT_A_dx, dS_dqU_B_bar );
          TensorDD term23 = outer( dT_A_dx_cen - dT_A_dx, tr_dS_bar ) / static_cast<double>(nDim);

          k_UU_mat += term23;
          k_UU_mat *= V0;

          // Geometrische Steifigkeit
          TensorDD k_UU_geo = ( - einsum< k, ij, i, to_jk >( dT_A_dx, S, dN_B_dx ) );
          TensorDD k_UU_geo_p = p_tau * outer(dN_B_dx, dT_A_dx) - p_tau * outer(dN_B_dx_cen, dT_A_dx_cen);

          k_UU_geo += k_UU_geo_p;
          k_UU_geo *= V0;

          TensorDD k_UU = k_UU_mat + k_UU_geo;
          k_UU += density0 * da_ddu * T_A * N_B * V0;

          {
              using namespace Eigen;
              K.template block< nDim, nDim >( idxA_u, idxB_u ) += Map< Matrix< double, nDim, nDim > >( torowmajor( k_UU ).data() );
          }
        }
      }
    }
  }
} // namespace Marmot::Meshfree
