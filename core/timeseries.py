import numpy as np

"""
timeseries.py

author: Jesse Berwald

opened: May 29, 2013

Simple class for manipulating time series data from neural network data.

(Replaces old, code-bloated version.)
"""

class Timeseries( object ):

    def __init__( self, data, sample_rate=1, timevec=None, interp=False ):
        """
        data : 1D array or list of data points

        Optional:
        ---------
        
        timevec : array or list of times. For use in interpolation
        (with interp==True) if they are unequally spaced.
        """
        self.data = np.asarray( data )
        self.sample_rate = sample_rate
        if self.sample_rate > 1:
            self.data = self.data[::self.sample_rate]

        if timevec is not None:
            self.tvec = timevec
        # if timevec and interp:
        #     #self.tvec = timevec
        #     pass ## do some cool interpolation like in the matlab 'bifurcation' code

    
            
    def _linear_spline( self ):
        """
        Compute a simple linear spline on the data points to create equally spaced data.

        %LINEAR_SPLINE     compute the piecewise linear interpolant associated 
        %                  with a given set of interpolating points and function 
        %                  values
        %
        %     calling sequences:
        %             ls = linear_spline ( xi, fi )
        %             linear_spline ( xi, fi )
        %
        %     inputs:
        %             xi      vector containing the interpolating points 
        %                     (must be sorted in ascending order, and equally spaced!)
        %             fi      vector containing function values
        %                     the i-th entry in this vector is the function
        %                     value associated with the i-th entry in the 'xi'
        %                     vector
        %
        %     output:
        %             ls      three column matrix containing the information which
        %                     defines the piecewise linear interpolant
        %                     - first column: interpolating points
        %                     - second column: function values
        %                     - third column: slopes of linear pieces
        %
        %     NOTE:
        %             to evaluate the piecewise linear interpolant apply the
        %             routine SPLINE_EVAL to the output of this routine 
        %   
        """
        nt = len( self.tvec )
        nx = len( self.data ) 

        assert nt == nx, "number of ordinates and number of function values must be equal!"

        diff = np.diff( self.data )
        slopes = np.diff( self.data ) / np.diff( self.tvec )

        return slopes
    
    def interpolator( self ):
        """
        Compute the linear interpolation of the function y = f(x) at the
        points x=self.tvec. x must be monotone (increasing or decreasing,
        doesn't matter).

        In-place computation that replaces the unevenly-spaced data
        and replaces both tvec and data with the interpolant.
        """
        assert len( self.data ) == len( self.tvec ), "x and y must be of the same length!"

        # compute the slopes between each (x_j,y_j) and (x_j+1,y_j+1)
        dy = self._linear_spline()

        print dy

        # subsample is necessary
        nx = len( self.tvec ) / self.sample_rate
        xmin = self.tvec.min()
        xmax = self.tvec.max()
        xi = np.linspace( xmin, xmax, nx );

        print "xi"
        print xi
        
        yi = np.zeros_like( xi )
    
        #% build lines across each (xi(k), xi(k+1)) interval
        for i in range( 1, nx ):
            idx = np.where( ( (xi >= self.tvec[i-1] ) & \
                            ( (xi <= self.tvec[i] ) ) ) )[0]
            # y = m(x - xi) + yi
            yi[ idx ] = dy[i-1] * ( xi[idx] - self.tvec[i-1] ) + self.data[i-1]
            # print "YI", idx,  yi[ idx ]

        self.data = np.asarray( yi )
        self.tvec = xi


if __name__ == "__main__":

    import matplotlib.pyplot as plt

    # Test the linear_interpolator.m functions
    xf = 10 * np.random.random( 50 ) #np.arange(0,10,0.05)
    xf.sort()
    yf = np.sin(2*np.pi*xf/5)+1;
    xp = np.arange(10)
    yp = np.sin(2*np.pi*xp/5)+1

    fig1 = plt.figure()
    ax1 = fig1.gca()
    ax1.plot( xf, yf, 'b-s', label='original, unevenly spaced' )
    ax1.plot( xp, yp, 'g--', label='original, evenly, sparsely spaced' )

    ts = Timeseries( yf, timevec=xf )
                    
    print 'Testing ts.interpolator()'
    ts.interpolator();

    print 'Plotting interpolated function...'
    ax1.plot( ts.tvec, ts.data, 'r-o', label='interpolated' )
    ax1.legend()
    plt.show()
    
