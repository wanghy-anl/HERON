# Copyright 2020, Battelle Energy Alliance, LLC
# ALL RIGHTS RESERVED

from .ExampleValidator import Example
from .FARMValidators import RefGov_SESBOP_W, Para_RefGov_SESBOPTES_MW

known = {
    'Example': Example,
    'RefGov_SESBOP_W': RefGov_SESBOP_W,
    'Para_RefGov_SESBOPTES_MW': Para_RefGov_SESBOPTES_MW,
    # ModelicaGoverner: TODO,
}

def get_class(typ):
  """
    Returns the requested dispatcher type.
    @ In, typ, str, name of one of the dispatchers
    @ Out, class, object, class object
  """
  return known.get(typ, None)
