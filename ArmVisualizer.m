classdef ArmVisualizer < matlab.System
    % Host-side visualization for Simulink
    
    properties
        N  = 3
        Nx = 25
        UseFlipZ = true
        ShowStateLines = true
    end

    properties(Access=private)
        Fig1
        Ax1
        SegLines
        SegTips
        Tip

        Fig2
        Ax2
        Cur
        Dots
        Lines
        Leg
        Inited = false
    end

    methods(Access=protected)
        function setupImpl(obj)
            % Figure 1
            obj.Fig1 = figure(1); clf(obj.Fig1);
            obj.Ax1 = axes(obj.Fig1); hold(obj.Ax1,'on'); grid(obj.Ax1,'on');
            axis(obj.Ax1,'equal');
            xlabel(obj.Ax1,'X'); ylabel(obj.Ax1,'Y'); zlabel(obj.Ax1,'Z');
            view(obj.Ax1,[11,13]);
            xlim(obj.Ax1,[-1 1]); ylim(obj.Ax1,[-1 1]); zlim(obj.Ax1,[-1.5 1.1]);

            obj.SegLines = gobjects(obj.N,1);
            obj.SegTips  = gobjects(obj.N,1);
            for s=1:obj.N
                obj.SegLines(s) = plot3(obj.Ax1, NaN,NaN,NaN, 'LineWidth',2);
                col = obj.SegLines(s).Color;
                obj.SegTips(s) = plot3(obj.Ax1, NaN,NaN,NaN, 's', ...
                    'MarkerFaceColor',col,'MarkerEdgeColor','none','MarkerSize',7);
            end
            obj.Tip = plot3(obj.Ax1, NaN,NaN,NaN,'o', ...
                'MarkerFaceColor',[.8 .2 .2],'MarkerEdgeColor','none','MarkerSize',7);

            % Figure 2
            obj.Fig2 = figure(2); clf(obj.Fig2);
            obj.Ax2 = axes(obj.Fig2); hold(obj.Ax2,'on'); grid(obj.Ax2,'on');
            xlabel(obj.Ax2,'time (s)'); ylabel(obj.Ax2,'length change (m)');
            title(obj.Ax2,sprintf('length change of %d Sections', obj.N));

            obj.Cur = plot(obj.Ax2, [0 0], [0 1], 'k--', 'LineWidth',1.2);

            obj.Dots = gobjects(2*obj.N,1);
            obj.Lines = gobjects(2*obj.N,1);

            % Pre-create dots; lines optional (lines require history buffering)
            for i=1:2*obj.N
                obj.Dots(i) = plot(obj.Ax2, NaN, NaN, 'o', ...
                    'MarkerFaceColor','auto','MarkerEdgeColor','none','MarkerSize',6);
            end

            % Legend labels
            leg = strings(2*obj.N,1);
            for ksec = 1:obj.N
                leg(2*(ksec-1)+1) = "l_{" + ksec + "1}";
                leg(2*(ksec-1)+2) = "l_{" + ksec + "2}";
            end
            legend(obj.Ax2, leg, 'Location','best');

            obj.Inited = true;
        end

        function stepImpl(obj, t, X, P)
            if ~obj.Inited || ~isvalid(obj.Fig1) || ~isvalid(obj.Fig2)
                setupImpl(obj);
            end

            % ---- Update arm (Figure 1) ----
            for s=1:obj.N
                Xs = P(:,1,s); Ys = P(:,2,s); Zs = P(:,3,s);
                set(obj.SegLines(s), 'XData', Xs, 'YData', Ys, 'ZData', Zs);

                % tip marker = last point
                set(obj.SegTips(s), 'XData', Xs(end), 'YData', Ys(end), 'ZData', Zs(end));
            end
            % overall tip = last section last point
            set(obj.Tip, 'XData', P(end,1,obj.N), 'YData', P(end,2,obj.N), 'ZData', P(end,3,obj.N));

            % ---- Update states (Figure 2) ----
            yl = ylim(obj.Ax2);
            set(obj.Cur, 'XData', [t t], 'YData', yl);

            for i=1:2*obj.N
                set(obj.Dots(i), 'XData', t, 'YData', X(i));
            end

            drawnow limitrate
        end

        function resetImpl(obj)
            % no-op
        end
    end
end
