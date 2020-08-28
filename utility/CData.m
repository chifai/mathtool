%% class to store and manipulate raw data
classdef CData < handle
    properties
        m_RawData = []
    end
    
    properties (SetAccess = private)
        m_Period = 0.01;    % default as 10 ms
        m_nDataSize;
        m_nDataType;
        m_TimeLine;
        m_cLegend;
        m_FigureNum = 0;
        m_xLabel = 'Time (s)'
        m_yLabel = 'Value'
        m_Title = 'Untitled'
        m_LineWidth = 1.5
    end
    
    methods
        function obj = CData(RawData, Period)
            obj.m_RawData = RawData;
            if nargin == 2
                obj.m_Period = Period;
            end
            % set properties
            obj.m_nDataSize = size(RawData, 1);
            obj.m_nDataType = size(RawData, 2);
            obj.m_TimeLine = 0:obj.m_Period:(obj.m_nDataSize-1) * obj.m_Period;
        end
        
        function Scale(obj, Factor)
            nFactorSize = size( Factor, 2 );
            for i = 1:obj.m_nDataType
                if i <= nFactorSize
                    obj.m_RawData(:, i) = obj.m_RawData(:, i) * Factor(i);
                else
                    break;
                end
            end
        end
        
        function SetLegend(obj, cLegend)
            obj.m_cLegend = cLegend;
        end
        
        function SetLabel(obj, xLabel, yLabel)
            obj.m_xLabel = xLabel;
            obj.m_yLabel = yLabel;
        end
        
        function SetTitle(obj, Title)
            obj.m_Title = Title;
        end
        
        function plot2D(obj)
            figure;
            obj.m_FigureNum = obj.m_FigureNum + 1;
            hold on;
            grid on;
            grid minor;
            for i = 1:obj.m_nDataType
                if strcmp(obj.m_cLegend{i}, '-') == 1
                    continue
                end
                plot(obj.m_TimeLine, obj.m_RawData(:,i), 'linewidth', obj.m_LineWidth);
                legend(obj.m_cLegend{i});
            end
            nValidIndex = find( strcmp(obj.m_cLegend, '-') == 0 );
            legend(obj.m_cLegend(nValidIndex));
            title(obj.m_Title);
            xlabel(obj.m_xLabel);
            ylabel(obj.m_yLabel);
            hold off;
        end
    end
end