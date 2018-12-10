classdef testGem < matlab.unittest.TestCase
    properties (TestParameter)
        type = {'single','double','uint16'};
        outSize = struct('s2d',[3 3], 's3d',[2 5 4]);
    end
    
    methods (Test)
        function testClass(testCase, type, outSize)
            testCase.verifyClass(zeros(outSize,type), type);
        end
        
        function testSize(testCase, outSize)
            testCase.verifySize(zeros(outSize), outSize);
        end
        
        function testDefaultClass(testCase)
            testCase.verifyClass(zeros, 'double');
        end
        function testDefaultSize(testCase)
            testCase.verifySize(zeros, [1 1]);
        end
        
        function testDefaultValue(testCase)
            testCase.verifyEqual(zeros,0);
        end
    end
end