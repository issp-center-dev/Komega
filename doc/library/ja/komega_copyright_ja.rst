プログラムの再配布
==================

自分のプログラムにKωを含める
-----------------------------

:math:`K\omega` ライブラリは下記の :ref:`lgplicense` (LGPL)に基づいて配布されている.
これはかいつまんで言うと次のようなことである.

 * 個人的なプログラムや, 研究室や共同研究者等のグループでソースコードをやりとりする時には,
   自由にコピペしたり改変して良い.
   
 * 公開したり売ったりするプログラムに関しては次のとおりである.
   
    * 配布するソースコードに :math:`K\omega` をそのまま,
      あるいは改変して含めるときには, そのプログラム本体をLGPL/GPLで配布する.
      
    * 配布するソースコードに含めず呼び出すだけならばライセンスによらず自由に配布できる.
      
    * ただしバイナリファイルを配布する場合に, そのバイナリに :math:`K\omega` が
      静的リンクされている場合にはLGPL/GPLで配布する.
      動的リンクされている(したがって :math:`K\omega` そのものはバイナリに含まれていない)
      場合にはライセンスによらず自由に配布できる.

.. _lgplicense
      
Lesser General Public License
-----------------------------

*© 2016- The University of Tokyo. All rights reserved.*

| This software is developed under the support of
| "*Project for advancement of software usability in materials science*" by The
| Institute for Solid State Physics, The University of Tokyo.
|
| This library is free software; you can redistribute it and/or
| modify it under the terms of the GNU Lesser General Public
| License as published by the Free Software Foundation; either
| version 2.1 of the License, or (at your option) any later version.
| This library is distributed in the hope that it will be useful,
| but WITHOUT ANY WARRANTY; without even the implied warranty of
| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
| Lesser General Public License for more details.
|
| You should have received a copy of the GNU Lesser General Public
| License along with this library; if not, write to the Free Software
| Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
|
| For more details, See 'COPYING.LESSER' in the root directory of this library.

