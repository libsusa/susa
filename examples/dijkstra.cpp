/*
 * This file is part of Susa.
 *
 * Susa is free software: you can redistribute it and/or modify
 * it under the terms of the Lesser GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * at your option) any later version.
 *
 * Susa is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * Lesser GNU General Public License for more details.
 *
 * You should have received a copy of the Lesser GNU General Public License
 * along with Susa.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <susa.h>
#include <limits> 
#include <algorithm>

const unsigned int num_nodes = 350;

int main(void)
{
  susa::matrix <unsigned int> graph(num_nodes, num_nodes, 0u);


  graph(10,5) = 12;graph(10,32) = 2;
  graph(5,21) = 4;graph(32,55) = 44;
  graph(21,55) = 3;


  susa::matrix <unsigned int> prev = dijkstra(graph, 10);

  unsigned int target = 55;
  unsigned int u = target;

  while (prev(u))
  {
    std::cout << " " << u;
    u = prev(u);
  }
  std::cout << " " << u << std::endl;

  return 0;
}
