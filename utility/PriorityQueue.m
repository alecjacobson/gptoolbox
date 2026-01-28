% There is some overhead using this class rather than just
% java.util.PriorityQueue and encode/decode directly.
classdef PriorityQueue < handle
  properties
    jQ
  end
  methods
    function pQ = PriorityQueue()
      pQ.jQ = java.util.PriorityQueue;
    end
    function pQ = add(pQ,value,index)
      pQ.jQ.add(PriorityQueue.encode(value,index));
    end
    function [value,index] = remove(pQ)
      [value,index] = PriorityQueue.decode(pQ.jQ.remove());
    end
    function res = isempty(pQ)
      res = pQ.jQ.isEmpty();
    end
    function [value,index] = peek(pQ)
      [value,index] = PriorityQueue.decode(pQ.jQ.peek());
    end
  end
  methods(Static)
    function u = encode(value,index)
      % value should be non-negative
      % index should be an integer
      u = sprintf('%019d,%d',typecast(value,'uint64'),index);
    end
    function [value,index] = decode(u)
      value_index = sscanf(u,'%lu,%lu',2);
      value = typecast(value_index(1),'double');
      index = double(value_index(2));
    end
  end
end
