classdef NoCTranspose < double
   methods
      function obj = ctranspose(data)
        error("You called ctranspose (') on a NoCTranspose object");
      end
   end
end
