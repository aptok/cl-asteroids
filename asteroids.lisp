(require 'magicl)
(require 'local-time)
(require 'trivial-download)

;; useful stuff
(defparameter *rad* (/ pi 180))
(defparameter *deg* (/ 180 pi))
(defun parse-float (string)
  "Return a float read from string, and the index to the remainder of string.
   Stolen from Peter Norvig"
  (multiple-value-bind (integer i)
      (parse-integer string :junk-allowed t)
    (multiple-value-bind (fraction j)
	(parse-integer string :start (+ i 1) :junk-allowed t)
      (values (float (+ integer (/ fraction (expt 10d0 (- j i 1))))) j))))

;; orbital elements of small bodies downloaded from
;; https://ssd.jpl.nasa.gov/dat/ELEMENTS.NUMBR and
;; https://ssd.jpl.nasa.gov/dat/ELEMENTS.UNNUM
(defparameter *elements-unnumbered-file* "ELEMENTS.UNNUM")
(defparameter *elements-numbered-file* "ELEMENTS.NUMBR")
(defparameter *accuracy* 1d-8) ; used in kepler-newton
;; astrodynamic parameters copied from
;; https://ssd.jpl.nasa.gov/astro_par.html
(defparameter *ecliptic-obliquity* (* *rad* (* (/ 1d0 3600d0) 84381.448d0))) ; radians
;(defparameter *J2000* 2451545d0) ; days since TODO
(defparameter *J2000* (local-time:encode-timestamp 0 0 0 1 1 1 2000)) ; J2000 timestamp
(defparameter *au* 149597870700) ; meters

;; Gaussian gravitational constant, not a constant anymore
;; https://www.iau.org/static/resolutions/IAU2012_English.pdf
(defparameter *k* 0.01720209895)

;; functions to read the orbital elements
(defun download-element-files-from-jpl ()
  "Download the orbital elements from JPL"
  (trivial-download:download "https://ssd.jpl.nasa.gov/dat/ELEMENTS.NUMBR" "ELEMENTS.NUMBR")
  (trivial-download:download "https://ssd.jpl.nasa.gov/dat/ELEMENTS.UNNUM" "ELEMENTS.UNNUM"))

(defun read-elements-format (format-string)
  "INPUT: second line of ELEMENTS.* file consisting of dashes and whitespaces.
   OUTPUT: positions of whitespaces in INPUT string"
  (let* ((p 0)
	 (b 0))
    (loop while (setf b (position #\  format-string :start (+ p 1)))
	  do (setf p b)
	     (setf b (position #\  format-string :start (+ p 1)))
		  collect p)))

(defun read-element-line-rec (format-list elements-line start ret)
  "INPUT: positions of white spaces as given by read-elements-format
          a line from ELEMENTS.* file
          start is for recursion
          ret is the accumulator
   OUTPUT: list of elements as strings form elements-line"
  (cond (format-list (setf ret (cons (subseq elements-line start (car format-list)) ret))
		     (read-element-line-rec (cdr format-list) elements-line (+ (car format-list) 1) ret))
	(t  (setf ret (cons (subseq elements-line start (car format-list)) ret)) (reverse ret))))

(defun read-element-line (format-list elements-line)
  "INPUT: positions of white spaces as given by read-elements-format
          a line from ELEMENTS.* file
   OUTPUT: list of elements as strings form elements-line"
  (read-element-line-rec format-list elements-line 0 ()))

;; functions for orbital mechanics
(defun kepler-newton (mean-anomaly excentricity excentric-anomaly-guessed accuracy)
  "INPUT: mean anomaly and a guessed excentric anomaly, given as radians
          excentricity
          desired accuracy of the result
   OUTPUT: excentric anomaly in radians
   Example: e=0.20563, M=3.3592, accuracy=5d-5, M as guess, should give 3.3222
   (kepler-newton 3.3592 0.20563 3.3592 5d-5)
   3.3222401615522115d0"
  (let ((excentric-anomaly-better (+ mean-anomaly (* excentricity (sin excentric-anomaly-guessed))))
	(excentricity (coerce excentricity 'double-float)))
    (cond ((>= (* (/ excentricity (- 1 excentricity)) (abs (- excentric-anomaly-better excentric-anomaly-guessed))) accuracy)
	   (kepler-newton mean-anomaly excentricity excentric-anomaly-better accuracy))
	  (t excentric-anomaly-better))))

(defun excentric-anomaly (excentricity mean-anomaly)
  "INPUT: excentricity, mean anomaly in radians
   OUTPUT: excentric anomaly in radians
   Example: e=0.20563 M=3.3592 should give
   (excentric-anomaly 0.20563 3.3592)
   3.32225263703146d0"
  (let ((m (mod mean-anomaly (* 2 pi))))
    (cond ((> m (* 2 pi))
	   (excentric-anomaly excentricity (- m (* 2 pi))))
	  ((< m (* -2 pi))
	   (excentric-anomaly excentricity (+ m (* 2 pi))))
	  (t 
	   (kepler-newton m excentricity m *accuracy*)))))

(defun cos-true-anomaly (excentricity excentric-anomaly)
  "INPUT: excentricity, excentric anomaly in radians
   OUTPUT: cosinus of the true anomaly"
  (let ((cos-excentric-anomaly (cos excentric-anomaly)))
    (/ (- cos-excentric-anomaly excentricity) (- 1 (* excentricity cos-excentric-anomaly)))))

(defun true-anomaly (e m)
  "INPUT: excentricity e, mean anomaly m
   OUTPUT: true anomaly
   Example: e=0.20563 M=9.6424 should give 2*pi-2.9948=3.2884
   (true-anomaly 0.20563 3.3592)
   3.2883734115413286d0"
  (let* ((excentric-anomaly (excentric-anomaly e m)))
    (cond ((and (<= 0 excentric-anomaly) (<= excentric-anomaly pi))
	   (acos (cos-true-anomaly e excentric-anomaly)))
	  ((and (<= pi excentric-anomaly) (<= excentric-anomaly (* 2 pi)))
	   (- (* 2 pi) (acos (cos-true-anomaly e excentric-anomaly))))
	  ((> excentric-anomaly (* 2 pi))
	   (error "excentric-anomaly bigger than 2*pi")))))


(defun distance-to-sun (a e m)
  "INPUT: semimajor axis a, excentricity e, mean anomaly m in radians
   OUTPUT: distance to sun in au
   Example: a=5.202921, e=0.04831, M=187°.5453=3.2732rad should give r=5.4523
   (distance-to-sun 5.202921 0.04831 3.2732)
   5.452295092212026d0"
  (/ (* a (- 1 (* e e)))
     (+ 1 (* e (cos-true-anomaly e (excentric-anomaly e m)))))) 

(defun d1 (alpha)
  "rotation matrix for the x axis"
  (magicl:from-list (list 1d0 0d0 0d0
			  0d0 (cos alpha) (- (sin alpha))
			  0d0 (sin alpha) (cos alpha))
		    '(3 3)))

(defun d2 (alpha)
  "rotation matrix for the y axis"
  (magicl:from-list (list (cos alpha) 0d0 (sin alpha)
			  0d0 1d0 0d0
			  (- (sin alpha)) 0d0 (cos alpha))
		    '(3 3)))

(defun d3 (alpha)
  "rotation matrix for the z axis"
  (magicl:from-list (list (cos alpha) (- (sin alpha)) 0d0
			  (sin alpha) (cos alpha) 0d0
			  0d0 0d0 1d0)
		    '(3 3)))

(defun s1 (phi)
  "vector in the orbital plane"
  (magicl:from-list (list (cos phi)
			  (sin phi)
			  0d0)
		    '(3 1)))

(defun s2 (phi)
  "vectro in the orbital plane"
  (magicl:from-list (list (- (sin phi))
			  (cos phi)
			  0d0)
		    '(3 1)))
  
(defun julian-centuries-from-j2000 (timestamp)
  "INPUT: timestamp as given by local-time package
   OUTPUT: julian centuries since/before 1. Jan. 2000"
  (/ (-
      (local-time:astronomical-julian-date timestamp)
      *J2000*)
     36525d0))

(defun degminsec2degdec (angle m s)
  "converts angle given in degree,minutes and seconds to decimal degrees"
  (+ angle (* (/ 1 60) m) (* (/ 1 3600) s)))

(defun transform-ecliptic-to-equatorial (pos)
  "INPUT: heliocentric ecliptic position
   OUTPUT: hecliocentric equatorial position"
  (magicl:@ (d1 *ecliptic-obliquity*) pos))

(defun mean-daily-movement (a)
  "INPUT: mean anomaly a in radians
   OUTPUT: mean daily movement in radians"
    (* *k* (sqrt (/ 1 (* a a a)))))

(defun heliocentric-position (a e m node i w)
  "INPUT: semimajor axis a, excentricity e, mena anomaly m,
          longitude of the ascending node node, inclination i,
          argument of periapsis w
          m, node, i and w in radians
   OUTPUT: hecliocentric ecliptic cordinates"
  (let* ((m (mod m (* 2 pi)))
	 (phi (true-anomaly e m))
	 (r0 (distance-to-sun a e m))
	 (d3node (d3 node))
	 (r0d3node (magicl:scale d3node r0))
	 (d1i (d1 i))
	 (d3o (d3 w))
	 (s1 (magicl:from-list (list (cos-true-anomaly e (excentric-anomaly e m)) (sin phi) 0) '(3 1))))
    (magicl:@ r0d3node d1i d3o s1)))

(defun heliocentric-position-at-date (a e m node i w epoch date)
    "INPUT: semimajor axis a, excentricity e, mena anomaly m,
          longitude of the ascending node node, inclination i,
          argument of periapsis w, modified julian date
          m, node, i and w in radians
   OUTPUT: hecliocentric ecliptic cordinates at specified date"
  (let* ((m0 (mod m (* 2 pi)))
	 (m1 (+ m0(* (mean-daily-movement a) (- date epoch)))))
	 (heliocentric-position a e m1 node i w)))

;; "class to store all the elements for an asteroid.
;;  unnumbered asteroids have identical name and id"
(defclass asteroid ()
  ((id :initarg :id :accessor id)
   (name :initarg :name :accessor name)
   (epoch :initarg :epoch :accessor epoch)
   (semimajor-axis :initarg :a :accessor a)
   (excentricity :initarg :e :accessor e)
   (inclination :initarg :i :accessor i)
   (argument-of-perihelion :initarg :w :accessor w)
   (longitude-of-ascending-node :initarg :node :accessor node)
   (mean-anomaly :initarg :m :accessor m)
   (absolute-magnitude :initarg :h :accessor h)
   (slope-parameter :initarg :g :accessor g)
   (ref :initarg :ref :accessor ref)))

(defun asteroid-mean-daily-movement (asteroid)
  "INPUT: asteroid object
   OUTPUT: mean daily movement in radians"
  (mean-daily-movement (a asteroid)))

(defun asteroid-heliocentric-position (asteroid)
   "INPUT: asteroid object
    OUTPUT: hecliocentric ecliptic cordinates"
  (heliocentric-position (a asteroid) (e asteroid)
			 (m asteroid) (node asteroid)
			 (i asteroid) (w asteroid)))

(defun asteroid-heliocentric-position-at-date (asteroid date)
     "INPUT: asteroid object, modified julian date
      OUTPUT: hecliocentric ecliptic cordinates"
  (heliocentric-position-at-date (a asteroid) (e asteroid)
				 (m asteroid) (node asteroid)
				 (i asteroid) (w asteroid)
				 (epoch asteroid) date))

(defun asteroid-from-numbered-elements-list (list)
  "INPUT: list of unnumbered asteroid orbital elements given by read element line
   OUTPUT: asteroid object"
  (let ((id (read-from-string (nth 0 list)))
	(name (string-trim '(#\Space #\Tab) (nth 1 list)))
	(epoch (read-from-string (nth 2 list)))
	(semi-major-axis (parse-float (nth 3 list)))
	(excentricity (parse-float (nth 4 list)))
	(inclination (* *rad* (parse-float (nth 5 list))))
	(argument-of-perihelion (* *rad* (parse-float (nth 6 list))))
	(longitude-of-ascending-node (* *rad* (parse-float (nth 7 list))))
	(mean-anomaly (* *rad* (parse-float (nth 8 list))))
	(absolute-magnitude (read-from-string (nth 9 list)))
	(slope-parameter (read-from-string (nth 10 list)))
	(ref (string-trim '(#\Space #\Tab) (nth 11 list))))
    (make-instance 'asteroid :id id :name name :epoch epoch
				      :a semi-major-axis :e excentricity 
				      :i inclination :w argument-of-perihelion
				      :node longitude-of-ascending-node :m mean-anomaly 
				      :h absolute-magnitude :g slope-parameter :ref ref)))

(defun asteroid-from-unnumbered-elements-list (list)
  "INPUT: list of numbered asteroid orbital elements given by read element line
   OUTPUT: asteroid object"
  (let* ((id (string-trim '(#\Space #\Tab) (nth 0 list)))
	 (name id)
	 (epoch (read-from-string (nth 1 list)))
	 (semi-major-axis (parse-float (nth 2 list)))
	 (excentricity (parse-float (nth 3 list)))
	 (inclination (* *rad* (parse-float (nth 4 list))))
	 (argument-of-perihelion (* *rad* (parse-float (nth 5 list))))
	 (longitude-of-ascending-node (* *rad* (parse-float (nth 6 list))))
	 (mean-anomaly (* *rad* (parse-float (nth 7 list))))
	 (absolute-magnitude (read-from-string (nth 8 list)))
	 (slope-parameter (read-from-string (nth 9 list)))
	 (ref (string-trim '(#\Space #\Tab) (nth 10 list))))
    (make-instance 'asteroid :id id :name name :epoch epoch
				      :a semi-major-axis :e excentricity 
				      :i inclination :w argument-of-perihelion
				      :node longitude-of-ascending-node :m mean-anomaly 
				      :h absolute-magnitude :g slope-parameter :ref ref)))

(defun position-key (pos)
  "INPUT: position vector as magicl object
   OUTPUT: hash key"
  (+ (floor (* 10 (+ 50 (magicl:tref pos 0 0))))
     (* 1000 (floor (* 10 (+ 50 (magicl:tref pos 1 0)))))
     (* 1000000 (floor (* 10 (+ 50 (magicl:tref pos 2 0)))))))

(defun asteroid-hash (asteroid hashtable)
  "INPUT: asteroid object, hash table
   OUTPUT: puts asteroid object into hash table"
  (let* ((pos (heliocentric-position (a asteroid) (e asteroid) (m asteroid)
				     (node asteroid) (i asteroid) (w asteroid)))
	 (key (position-key pos)))
    (setf (gethash key hashtable) (cons pos (gethash key hashtable)))))
	
(defun read-elements-numbered ()
  "INPUT: file with numbered elements from https://ssd.jpl.nasa.gov/dat/ELEMENTS.NUMBR
   OUTPUT: list of asteroid elements"
  (with-open-file (in *elements-numbered-file*)
    (read-line in)
    (let* ((format-string (read-line in))
	   (format-list (read-elements-format format-string)))
      (loop for line = (read-line in nil)
            while line
            collect (asteroid-from-numbered-elements-list (read-element-line format-list line))))))

(defun read-elements-unnumbered ()
    "INPUT: file with numbered elements from https://ssd.jpl.nasa.gov/dat/ELEMENTS.UNNUM
     OUTPUT: list of asteroid objects"
  (with-open-file (in *elements-unnumbered-file*)
    (read-line in)
    (let* ((format-string (read-line in))
	   (format-list (read-elements-format format-string)))
      (loop for line = (read-line in nil)
            while line
            collect (asteroid-from-unnumbered-elements-list (read-element-line format-list line))))))

(defun hash-elements (hashtable asteroid-list)
  "INPUT: hash table, list of asteroid elementss from read-elements-numbered/unnumbered
   OUTPUT: fills hashtable with the asteroid objects"
  (dolist (asteroid asteroid-list)
    (asteroid-hash asteroid hashtable)))

(defun unpack-magicl-position (pos)
  (list (magicl:tref pos 0 0) (magicl:tref pos 1 0) (magicl:tref pos 2 0)))
