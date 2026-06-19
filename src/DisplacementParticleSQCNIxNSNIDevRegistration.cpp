#include "Marmot/DisplacementParticleSQCNIxNSNIDev.h"
#include "Marmot/MarmotParticleLibrary.h"

namespace Marmot::Meshfree {
  using namespace MarmotLibrary;

  // Registration for PlaneStrain/Quad (2D, 4 vertices)
  const static bool DisplacementParticleSQCNIxNSNIDev_PlaneStrain_Quad_isRegistered =
    MarmotLibrary::MarmotParticleFactory::registerParticle(
      "DisplacementSQCNIxNSNIDev/PlaneStrain/Quad",
      []( int cellID,
          const double* nodeCoordinates,
          int sizeNodeCoordinates,
          double volume,
          const std::string& materialName,
          const double* materialProperties,
          int sizeMaterialProperties,
          const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation ) -> Marmot::Meshfree::MarmotParticle* {
        return new DisplacementParticleSQCNIxNSNIDev< 2, 4 >(
          cellID,
          nodeCoordinates,
          sizeNodeCoordinates,
          volume,
          materialName,
          materialProperties,
          sizeMaterialProperties,
          approximation,
          DisplacementParticleSQCNIxNSNIDev< 2, 4 >::SmoothingDomainUpdateType::DeformationGradient );
      } );

  // Registration for 3D/Hexa (3D, 8 vertices)
  const static bool DisplacementParticleSQCNIxNSNIDev_3D_Hexa_isRegistered =
    MarmotLibrary::MarmotParticleFactory::registerParticle(
      "DisplacementSQCNIxNSNIDev/3D/Hexa",
      []( int cellID,
          const double* nodeCoordinates,
          int sizeNodeCoordinates,
          double volume,
          const std::string& materialName,
          const double* materialProperties,
          int sizeMaterialProperties,
          const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation ) -> Marmot::Meshfree::MarmotParticle* {
        return new DisplacementParticleSQCNIxNSNIDev< 3, 8 >(
          cellID,
          nodeCoordinates,
          sizeNodeCoordinates,
          volume,
          materialName,
          materialProperties,
          sizeMaterialProperties,
          approximation,
          DisplacementParticleSQCNIxNSNIDev< 3, 8 >::SmoothingDomainUpdateType::DeformationGradient );
      } );

  // Registration for PlaneStrain/Quad (2D, 4 vertices) with Nonconforming update
  const static bool DisplacementParticleSNNIxNSNIDev_PlaneStrain_Quad_isRegistered =
    MarmotLibrary::MarmotParticleFactory::registerParticle(
      "DisplacementSNNIxNSNIDev/PlaneStrain/Quad",
      []( int cellID,
          const double* nodeCoordinates,
          int sizeNodeCoordinates,
          double volume,
          const std::string& materialName,
          const double* materialProperties,
          int sizeMaterialProperties,
          const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation ) -> Marmot::Meshfree::MarmotParticle* {
        return new DisplacementParticleSQCNIxNSNIDev< 2, 4 >(
          cellID,
          nodeCoordinates,
          sizeNodeCoordinates,
          volume,
          materialName,
          materialProperties,
          sizeMaterialProperties,
          approximation,
          DisplacementParticleSQCNIxNSNIDev< 2, 4 >::SmoothingDomainUpdateType::None);
      } );

  // Registration for 3D/Hexa (3D, 8 vertices) with nonconforming update
  const static bool DisplacementParticleSNNIxNSNIDev_3D_Hexa_isRegistered =
    MarmotLibrary::MarmotParticleFactory::registerParticle(
      "DisplacementSQCNIxNSNIDev/3D/Hexa",
      []( int cellID,
          const double* nodeCoordinates,
          int sizeNodeCoordinates,
          double volume,
          const std::string& materialName,
          const double* materialProperties,
          int sizeMaterialProperties,
          const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation ) -> Marmot::Meshfree::MarmotParticle* {
        return new DisplacementParticleSQCNIxNSNIDev< 3, 8 >(
          cellID,
          nodeCoordinates,
          sizeNodeCoordinates,
          volume,
          materialName,
          materialProperties,
          sizeMaterialProperties,
          approximation,
          DisplacementParticleSQCNIxNSNIDev< 3, 8 >::SmoothingDomainUpdateType::None);
      } );

  const static bool
    DisplacementParticleSQCNIxNSNIDev_R_PlaneStrain_Quad_isRegistered = MarmotLibrary::MarmotParticleFactory::
      registerParticle( "DisplacementSQCNI_RxNSNIDev/PlaneStrain/Quad",
                        []( int                                                  cellID,
                            const double*                                        nodeCoordinates,
                            int                                                  sizeNodeCoordinates,
                            double                                               volume,
                            const std::string&                                   materialName,
                            const double*                                        materialProperties,
                            int                                                  sizeMaterialProperties,
                            const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                          -> Marmot::Meshfree::MarmotParticle* {
                          return new DisplacementParticleSQCNIxNSNIDev<
                            2,
                            4 >( cellID,
                                 nodeCoordinates,
                                 sizeNodeCoordinates,
                                 volume,
                                 materialName,
                                 materialProperties,
                                 sizeMaterialProperties,
                                 approximation,
                                 DisplacementParticleSQCNIxNSNIDev< 2, 4 >::SmoothingDomainUpdateType::RotationOnly );
                        } );

  const static bool DisplacementParticleSQCNIxNSNIDev_RU_PlaneStrain_Quad_isRegistered = MarmotLibrary::
    MarmotParticleFactory::registerParticle( "DisplacementSQCNI_RUxNSNIDev/PlaneStrain/Quad",
                                             []( int                cellID,
                                                 const double*      nodeCoordinates,
                                                 int                sizeNodeCoordinates,
                                                 double             volume,
                                                 const std::string& materialName,
                                                 const double*      materialProperties,
                                                 int                sizeMaterialProperties,
                                                 const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                                               -> Marmot::Meshfree::MarmotParticle* {
                                               return new DisplacementParticleSQCNIxNSNIDev<
                                                 2,
                                                 4 >( cellID,
                                                      nodeCoordinates,
                                                      sizeNodeCoordinates,
                                                      volume,
                                                      materialName,
                                                      materialProperties,
                                                      sizeMaterialProperties,
                                                      approximation,
                                                      DisplacementParticleSQCNIxNSNIDev< 2, 4 >::
                                                        SmoothingDomainUpdateType::RotationAndPrincipalStretch );
                                             } );
} // namespace Marmot::Meshfree
