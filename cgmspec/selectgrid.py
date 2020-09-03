import numpy as np
import pylab as plt


def get_cells(h,incli,D,alpha,size):
    m = -np.tan(np.radians(incli))
    y0 = D*np.sin(np.radians(alpha))/np.cos(np.radians(incli))
    n = -m*y0
    y1 = ((h/2)-n)/m
    y2 = (-(h/2)-n)/m

    z1 = h/2
    z2 = -h/2
    print('y1,y2', y1, y2)
    t = 0
    y = y1
    z = z1
    loslenght = np.sqrt(((y2-y1)**2)+(z2-z1)**2)
    if y1<0:
        nygrill = int(y1/size)
    else:
        nygrill = int((y1/size)-1)

    nzgrill = int((h/2)/size)
    print('grillas', nygrill, nzgrill)

    cellnum = 0
    cellpos = []
    fig = plt.figure()
    ax = fig.gca()
    major_ticks = np.arange(-50, 50, size)
    ylos = np.linspace(y1,y2, 100)
    zlos = n + m*ylos
    #ax.set_yticks(np.arange(z1,z2,size))
    ax.set_yticks(major_ticks)
    ax.set_xticks(major_ticks)
    ax.grid(which='both')
    plt.plot(ylos,zlos)



    while t<loslenght:
          print('grillas', nygrill, nzgrill)
          ypos = (size*nygrill)+(size/2)
          zpos = (size*nzgrill)-(size/2)
          print('posicion', ypos,zpos)
          cellpos.append([ypos,zpos])
          cellnum = cellnum+1
          nexty = ((nygrill+1)*size, n+(m*(nygrill+1)*size))
          dnexty = np.sqrt(((nexty[0]-y)**2)+ (nexty[1]-z)**2)
          nextz = ((((nzgrill-1)*size)-n)/m, (nzgrill-1)*size)
          dnextz = np.sqrt(((nextz[0]-y)**2)+ (nextz[1]-z)**2)
          #nextycross = (size*nycell, n+(m*size*nycell))
          #nextzcross = ((z-n) / m, size*nzcell)

          #zdisttoycross = n+(m*y)
          #ydisttozcross = (z-n) / m
          #disttoycross = np.sqrt(((nextycross[0]-y)**2)+(nextycross[1]-z)**2)
          #disttozcross = np.sqrt(((nextzcross[0]-y)**2)+(nextzcross[1]-z)**2)
          #print('distancias', disttoycross, disttozcross)

          if dnexty < dnextz:
              print(0)
              t = t + dnexty
              y = nexty[0]
              z = nexty[1]
              nygrill = nygrill+1

          else:
              print(1)
              t = t + dnextz
              y = nextz[0]
              z = nextz[1]
              nzgrill = nzgrill-1


    for i in range(len(cellpos)):
        plt.plot(cellpos[i][0], cellpos[i][1], 'or')
    return(cellnum, cellpos)
