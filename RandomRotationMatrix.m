function [ R ] = RandomRotationMatrix()
    
    % Each column of a rotation matrix can be thought of as the x, y and
    % z-axis of some reference frame placed with some orientation.
    
    % We set x-axis to be a random 3D unit vector.
    % For it to have a "random" direction, we generate a random point
    % inside a unit sphere and use its direction.
    
    while (true)
        p = 2*(rand(3,1)-0.5);
        if (sum(p.*p) < 1)
            break;
        end
    end
    x_axis = p/norm(p);
    
    % Then, we generate y-axis by finding a vector that is perpendicular to
    % the x-axis. This is done by generating another random 3D unit vector
    % and then computing the cross product with the x-axis.
    
    while (true)
        p = 2*(rand(3,1)-0.5);
        if (sum(p.*p) < 1)
            break;
        end
    end
    
    y_axis = cross(x_axis, p);
    y_axis = y_axis/norm(y_axis);
    
    % z-axis must be the cross product of x-axis and y-axis.
    z_axis = cross(x_axis, y_axis);
    z_axis = z_axis/norm(z_axis);
    
    R = [x_axis, y_axis, z_axis];


end

