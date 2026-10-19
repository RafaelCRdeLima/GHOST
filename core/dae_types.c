#include "dae_types.h"

const char *dae_strerror(dae_status s)
{
  switch (s) {
    case DAE_OK:                 return "ok";
    case DAE_ERR_ALLOC:          return "memoria insuficiente";
    case DAE_ERR_PARAM:          return "argumento invalido";
    case DAE_ERR_JSON:           return "spec.json malformado";
    case DAE_ERR_VERSION:        return "versao de formato incompativel";
    case DAE_ERR_TOOBIG:         return "pedido acima do teto";
    case DAE_ERR_CANCELLED:      return "cancelado";
    case DAE_ERR_NORM:           return "norma nao conservada: limite espectral insuficiente";
    case DAE_ERR_UNIMPLEMENTED:  return "nao implementado";
  }
  return "erro desconhecido";
}
