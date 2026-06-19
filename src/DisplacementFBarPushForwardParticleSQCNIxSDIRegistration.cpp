#include "Marmot/DisplacementFBarPushForwardParticleSQCNIxSDI.h"
#include "Marmot/MarmotMaterialPoint.h"
#include "Marmot/MarmotParticleLibrary.h"

namespace Marmot::Meshfree {

  using namespace MarmotLibrary;

  const static bool DisplacementFBarPushForwardParticleSQCNIxSDI_PlaneStrain_Quad_isRegistered = MarmotLibrary::MarmotParticleFactory::
    registerParticle( "DisplacementFBarPushForwardSQCNIxSDI/PlaneStrain/Quad",
                      []( int                                                  cellID,
                          const double*                                        nodeCoordinates,
                          int                                                  sizeNodeCoordinates,
                          double                                               volume,
                          const std::string&                                   materialName,
                          const double*                                        materialProperties,
                          int                                                  sizeMaterialProperties,
                          const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                        -> Marmot::Meshfree::MarmotParticle* {
                        return new DisplacementFBarPushForwardParticleSQCNIxSDI<
                          2,
                          4 >( cellID,
                               nodeCoordinates,
                               sizeNodeCoordinates,
                               volume,
                               materialName,
                               materialProperties,
                               sizeMaterialProperties,
                               approximation,
                               DisplacementFBarPushForwardParticleSQCNIxSDI< 2, 4 >::SmoothingDomainUpdateType::DeformationGradient );
                      } );

  const static bool DisplacementFBarPushForwardParticleSQCNIxSDI_3D_Hexa_isRegistered = MarmotLibrary::MarmotParticleFactory::
    registerParticle( "DisplacementFBarPushForwardSQCNIxSDI/3D/Hexa",
                      []( int                                                  cellID,
                          const double*                                        nodeCoordinates,
                          int                                                  sizeNodeCoordinates,
                          double                                               volume,
                          const std::string&                                   materialName,
                          const double*                                        materialProperties,
                          int                                                  sizeMaterialProperties,
                          const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                        -> Marmot::Meshfree::MarmotParticle* {
                        return new DisplacementFBarPushForwardParticleSQCNIxSDI<
                          3,
                          8 >( cellID,
                               nodeCoordinates,
                               sizeNodeCoordinates,
                               volume,
                               materialName,
                               materialProperties,
                               sizeMaterialProperties,
                               approximation,
                               DisplacementFBarPushForwardParticleSQCNIxSDI< 3, 8 >::SmoothingDomainUpdateType::DeformationGradient );
                      } );

  const static bool DisplacementFBarPushForwardParticleRS_SNNIxSDI_3D_Hexa_isRegistered = MarmotLibrary::MarmotParticleFactory::
    registerParticle( "DisplacementFBarPushForwardRS-SNNIxSDI/3D/Hexa",
                      []( int                                                  cellID,
                          const double*                                        nodeCoordinates,
                          int                                                  sizeNodeCoordinates,
                          double                                               volume,
                          const std::string&                                   materialName,
                          const double*                                        materialProperties,
                          int                                                  sizeMaterialProperties,
                          const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                        -> Marmot::Meshfree::MarmotParticle* {
                        return new DisplacementFBarPushForwardParticleSQCNIxSDI< 3, 8 >( cellID,
                                                                          nodeCoordinates,
                                                                          sizeNodeCoordinates,
                                                                          volume,
                                                                          materialName,
                                                                          materialProperties,
                                                                          sizeMaterialProperties,
                                                                          approximation,
                                                                          DisplacementFBarPushForwardParticleSQCNIxSDI< 3, 8 >::
                                                                            SmoothingDomainUpdateType::
                                                                              RotationAndPrincipalStretch );
                      } );

  const static bool DisplacementFBarPushForwardParticleR_SNNIxSDI_3D_Hexa_isRegistered = MarmotLibrary::MarmotParticleFactory::
    registerParticle( "DisplacementFBarPushForwardR-SNNIxSDI/3D/Hexa",
                      []( int                                                  cellID,
                          const double*                                        nodeCoordinates,
                          int                                                  sizeNodeCoordinates,
                          double                                               volume,
                          const std::string&                                   materialName,
                          const double*                                        materialProperties,
                          int                                                  sizeMaterialProperties,
                          const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                        -> Marmot::Meshfree::MarmotParticle* {
                        return new DisplacementFBarPushForwardParticleSQCNIxSDI< 3, 8 >( cellID,
                                                                          nodeCoordinates,
                                                                          sizeNodeCoordinates,
                                                                          volume,
                                                                          materialName,
                                                                          materialProperties,
                                                                          sizeMaterialProperties,
                                                                          approximation,
                                                                          DisplacementFBarPushForwardParticleSQCNIxSDI< 3, 8 >::
                                                                            SmoothingDomainUpdateType::RotationOnly );
                      } );

  const static bool DisplacementFBarPushForwardParticleSNNIxSDI_3D_Hexa_isRegistered = MarmotLibrary::MarmotParticleFactory::
    registerParticle( "DisplacementFBarPushForwardSNNIxSDI/3D/Hexa",
                      []( int                                                  cellID,
                          const double*                                        nodeCoordinates,
                          int                                                  sizeNodeCoordinates,
                          double                                               volume,
                          const std::string&                                   materialName,
                          const double*                                        materialProperties,
                          int                                                  sizeMaterialProperties,
                          const Marmot::Meshfree::MarmotMeshfreeApproximation& approximation )
                        -> Marmot::Meshfree::MarmotParticle* {
                        return new DisplacementFBarPushForwardParticleSQCNIxSDI< 3, 8 >( cellID,
                                                                          nodeCoordinates,
                                                                          sizeNodeCoordinates,
                                                                          volume,
                                                                          materialName,
                                                                          materialProperties,
                                                                          sizeMaterialProperties,
                                                                          approximation,
                                                                          DisplacementFBarPushForwardParticleSQCNIxSDI< 3, 8 >::
                                                                            SmoothingDomainUpdateType::None );
                      } );

} // namespace Marmot::Meshfree
