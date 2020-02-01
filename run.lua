#!/usr/bin/env luajit
local table = require 'ext.table'
local class = require 'ext.class'
local math = require 'ext.math'
local View = require 'glapp.view'
local Orbit = require 'glapp.orbit'
local ImGuiApp = require 'imguiapp'
local matrix = require 'matrix'
local ffi = require 'ffi'
local vec3d = require 'vec-ffi.vec3d'
local ig = require 'ffi.imgui'
local sdl = require 'ffi.sdl'
local gl = require 'gl'

local App = class(Orbit(View.apply(ImGuiApp)))
App.title = 'wave in curved space'

-- collections of pts to render a lightcone
local rings

local pts
local n = 10

ffi.cdef[[
typedef struct {
	vec3d_t pos;
	vec3d_t vel;
} pt_t;
]]

local pt_t
pt_t = ffi.metatype('pt_t', {
	__add = function(a,b)
		return pt_t(
			a.pos + b.pos,
			a.vel + b.vel
		)
	end,
	__mul = function(a,s)
		return pt_t(
			a.pos * s,
			a.vel * s
		)
	end,
	__tostring = function(a)
		return tolua{pos=a.pos, vel=a.vel}
	end,
})

local minAngle = 0
local maxAngle = 0
local minDist = 0
local maxDist = 0

local body = {
	centerX = 0,
	centerY = 5,
	R = 1,	-- Schwarzschild radius
	rs = 1,	-- physical radius
}

function App:reset()
	self.t = 0
	pts = table()
	local v = 1
	local gamma = 1	-- / math.sqrt(1 - v * v)
	for i=0,n-1 do
		local theta = 2 * math.pi * i / n
		pts:insert(pt_t(
			vec3d(self.t, 0,0),
			vec3d(1, v * math.cos(theta), v * math.sin(theta)) * gamma
		))
	end

	rings = table()
end

App.updateMethod = false
App.t = 0
function App:init(...)
	App.super.init(self, ...)

	--self.view.ortho = true

	self:reset()
end

local function Conn_Vacuum(pt)
	return matrix.zeros(3,3,3)
end

local function Conn_Newton(pt)
	local d = pt - vec3d(0, body.centerX, body.centerY)
	local t,x,y = d:unpack()
	local r2 = x*x + y*y
	local epsilon = 2 * body.R
	local Conn = matrix.zeros(3,3,3)
	if r2 < epsilon * epsilon then return Conn end
	local r = math.sqrt(r2)
	local r_rs = r / body.rs
	local R = body.R * math.ceil(1, r_rs * r_rs)
	local r3 = r2 * r
	Conn[2][1][1] = R * x / r3
	Conn[3][1][1] = R * y / r3
	return Conn
end

--[=[
g_tt = -((1-R/(4r))/(1+R/(4r)))^2
g_ij = delta_ij (1+R/(4r))^4

g^tt = -((1+R/(4r))/(1-R/(4r)))^2
g^ij = delta^ij (1+R/(4r))^(-4)

g_tt,k = -2 ((1-R/(4r))/(1+R/(4r))) * (
			(1-R/(4r))' * (1+R/(4r)) - (1-R/(4r)) (1+R/(4r))'
		) / (1+R/(4r))^2
	= -(1-R/(4r))/(1+R/(4r))^3 * R/r^3 x^k

g_ij,k = delta_ij 4 (1+R/(4r))^3 (1+R/(4r))'
	= -R/r^3 (1+R/(4r))^3 delta_ij x^k

Conn_ttt = 0
Conn_itt = 1/2 (g_it,t + g_it,t - g_tt,i) = -1/2 g_tt,i = 1/2 R/r^3 (1-R/(4r))/(1+R/(4r))^3 x^l delta_li
Conn_tit = Conn_tti = 1/2 g_tt,i = -Conn_itt
Conn_ijt = Conn_itj = 1/2 (g_ij,t + g_it,j - g_jt,i) = 1/2 g_ij,t = 0 
Conn_tij = 1/2 (g_ti,j + g_tj,i - g_ij,t) = -1/2 g_ij,t = 0
Conn_ijk = 1/2 (g_ij,k + g_ik,j - g_jk,i) = 1/2 (
		- R/r^3 (1+R/(4r))^3 delta_ij x^k
		- R/r^3 (1+R/(4r))^3 delta_ik x^j
		+ R/r^3 (1+R/(4r))^3 delta_jk x^i
	)
	= R/(2 (r + R/4)^3) (-delta_ij delta_lk x^l - delta_ik delta_lj x^l + delta_jk delta_li x^l)

Conn^t_tt = g^tt Conn_ttt + g^tm Conn_mtt = 0
Conn^i_tt = g^it Conn_ttt + g^im Conn_mtt = delta^im (1+R/(4r))^(-4) * 1/2 R/r^3 (1-R/(4r))/(1+R/(4r))^3 x^l delta_lm
	= 1/2 R/r^3 (1-R/(4r))/(1+R/(4r))^7 x^i

Conn^t_it = Conn^t_ti = g^tt Conn_tti + g^tm Conn_mit = ((1+R/(4r))/(1-R/(4r)))^2 1/2 R/r^3 (1-R/(4r))/(1+R/(4r))^3 x^l delta_li
	= 1/2 R/r^3 1/( (1-R/(4r)) (1+R/(4r)) ) x^l delta_li

Conn^i_jt = Conn^i_tj = g^it Conn_ttj + g^im Conn_mjt = 0

Conn^t_ij = g^tt Conn_tij + g^tm Conn_mij = 0

Conn^i_jk = g^it Conn_tjk + g^im Conn_mjk = delta^im (1+R/(4r))^(-4) * 1/2 R/r^3 (1+R/(4r))^3 (-delta_mj delta_lk x^l - delta_mk delta_lj x^l + delta_jk delta_lm x^l)
	= 1/2 R/r^3 1/(1+R/(4r)) (-delta^i_j delta_lk x^l - delta^i_k delta_lj x^l + delta_jk x^i)

--]=]
local function Conn_Schwarzschild(pt)
	local d = pt - vec3d(0, body.centerX, body.centerY)
	local t,x,y = d:unpack()
	local r2 = x*x + y*y
	local epsilon = 2 * body.R
	local Conn = matrix.zeros(3,3,3)
	if r2 < epsilon * epsilon then return Conn end
	local _1_r2 = 1 / r2
	local r = math.sqrt(r2)
	
	local r_rs = r / body.rs
	local R = body.R * math.ceil(1, r_rs * r_rs)
	local R_r = R / r
	local R_r3 = R * _1_r2

	local mu_plus = 1 + .25 * R_r
	local mu_minus = 1 - .25 * R_r

	local mu_plus_2 = mu_plus * mu_plus   
	local mu_plus_4 = mu_plus_2 * mu_plus_2
	local mu_plus_7 = mu_plus_4 * mu_plus_2 * mu_plus 

	local conn_tti = .5 * R_r3 * 1 / (mu_plus * mu_minus)
	local conn_ttx = conn_tti * x
	Conn[1][1][2] = conn_ttx
	Conn[1][2][1] = conn_ttx
	local conn_tty = conn_tti * y 
	Conn[1][1][3] = conn_tty
	Conn[1][3][1] = conn_tty
	
	local conn_itt = .5 * R_r3 * mu_minus / mu_plus_7
	Conn[2][1][1] = conn_itt * x
	Conn[3][1][1] = conn_itt * y
	
	for i=1,2 do
		for j=1,2 do
			for k=1,2 do
				local s = 0
				if i == j then s = s - d.s[k] end
				if i == k then s = s - d.s[j] end
				if j == k then s = s + d.s[i] end
				Conn[i+1][j+1][k+1] = .5 * R_r3 / mu_plus * s
			end
		end
	end
	return Conn
end

--[[
g_uv = eta_uv + f k_u k_v
f = r^2 / (r^4 + a^2 z^2) (R r - Q^2)
k_u = (1, (rx+ay)/(r^2 + a^2), (ry-ax)/(r^2 + a^2), z/r)
a = J / M

rc^2 = x^2 + y^2 + z^2
1 = (x^2 + y^2) / (r^2 + a^2) + z^2 / r^2
r^2 (r^2 + a^2) - r^2 (x^2 + y^2) - z^2 (r^2 + a^2) = 0
r^4 + r^2 (a^2 - rc^2) - a^2 z^2 = 0
r = +- sqrt (
	rc^2 - a^2 +- sqrt(
		(rc^2 - a^2)^2 + 4 a^2 z^2
	)
)

... in 2D z=0:

f = R / r - Q^2 / r^2
f,u = -(R + 2 Q^2 / r) r,u / r^2

k_u = (1, (rx+ay)/(r^2 + a^2), (ry-ax)/(r^2 + a^2))
k_t,u = 0
k_x,x = ((r,x x + r) (r^2 + a^2) - (rx+ay) (2 r r,x)) / (r^2 + a^2)^2
k_y,x = ((r,x y - a) (r^2 + a^2) - (ry-ax) (2 r r,x)) / (r^2 + a^2)^2
k_x,y = ((r,y x + a) (r^2 + a^2) - (rx+ay) (2 r r,y)) / (r^2 + a^2)^2
k_y,y = ((r,y y + r) (r^2 + a^2) - (ry-ax) (2 r r,y)) / (r^2 + a^2)^2

rc^2 = x^2 + y^2
rc,a = x^a / sqrt(x^2 + y^2)
rc,a = x^a / rc

r = +-sqrt( rc^2 - a^2 +- |rc^2 - a^2| )
r = 0, r = +-sqrt(2 (rc^2 - a^2) )
r,a = +-(2 rc rc,a) / sqrt(2 (rc^2 - a^2) )
r,a = +- 2 x^a / sqrt(2(rc^2 - a^2))

g_ab,c = k_a k_b f,c + f (k_a,c k_b + k_a k_b,c)

--]]
local function Conn_KerrSchild(pt)
	local d = pt - vec3d(0, body.centerX, body.centerY)
	local t,x,y = d:unpack()
	local r2 = x*x + y*y
	local epsilon = 2 * body.R
	local Conn = matrix.zeros(3,3,3)
	if r2 < epsilon * epsilon then return Conn end
	local r = math.sqrt(r2)
	local rs = body.rs
	local r_rs = r / rs
	local R = body.R * math.ceil(1, r_rs * r_rs)

	local Q = body.Q or 0
	local a = body.a or 0	-- angular momentum (per unit mass?)


	-- r^4 + r^2 (a^2 - x^2 - y^2 - z^2) - a^2 z^2 = 0
	local rcSq = x*x + y*y
	local aSq = a*a
	local rSq = math.sqrt(rcSq - aSq)
	local r = math.sqrt(rSq)
	local H = .5 * (r * rs - Q * Q) / (rSq + aSq*z*z/rSq)
	local lU = matrix{
		-1,
		(r*x + a*y)/(rSq + a*a),
		(r*y - a*x)/(rSq + a*a),
		--z/r,
	}
	local gU = matrix{{-1,0,0},{0,1,0},{0,0,1}} - 2 * H * matrix.outer(lU,lU)


	local dr_dxis = matrix()
	for i=1,2 do
		dr_dxis[i] = (
			-.25 * db_dxis[i] 
			+ .5 * (
				b * db_dxis[i] 
				+ (i == 3 and (4*aSq*z) or 0)
			) / sqrtdiscr
		) / r
	end
	local dH_dx = matrix()
	dH_dx[1] = 0
	for i=1,2 do
		dH_dx[i+1] = .5 * (
			dr_dxis[i] * rs * (rSq + aSq*z*z / rSq)
			- (r * rs - Q*Q) * (
				2 * r * dr_dxis[i]
				+ aSq*(
					(i == 3 and (2 * z / rSq) or 0)
					- 2 * z*z*r*dr_dxis[i]
				)
			)
		)
	end

	local aSq_plus_rSq = rSq + aSq
	
	local l = matrix{
		1,
		(r*x + a*y)/aSq_plus_rSq,
		(r*y - a*x)/aSq_plus_rSq,
	}

	-- l_a,b == dl_dx[a][b]
	local dl_dx = matrix{
		{0,0,0},
		{
			0, 
			(dr_dxis[1] * (x * (a^2 - r^2) - 2 * a * y * r) / aSq_plus_rSq + r) / aSq_plus_rSq,
			(dr_dxis[2]  * (x * (a^2 - r^2) - 2 * a * y * r) / aSq_plus_rSq + a) / aSq_plus_rSq,
		},
		{
			0,
			(dr_dxis[1] * (y * (a^2 - r^2) + 2 * a * x * r) / aSq_plus_rSq - a) / aSq_plus_rSq,
			(dr_dxis[2] * (y * (a^2 - r^2) + 2 * a * x * r) / aSq_plus_rSq + r) / aSq_plus_rSq,
		},
	}

	-- conn_abc = H,a l_b l_c + H,b l_a l_c + H,c l_a l_b
	--		+ H (l_a (l_b,c + l_c,b) + l_b (l_a,c + l_c,a) + l_c (l_a,b + l_b,a))
	local connL = matrix{3,3,3}:lambda(function(a,b,c)
		return dH_dx[a] * l[b] * l[c]
			+ dH_dx[b] * l[a] * l[c]
			+ dH_dx[c] * l[a] * l[b]
			+ H * (
				l[a] * (dl_dx[b][c] + dl_dx[c][b])
				+ l[b] * (dl_dx[a][c] + dl_dx[c][a])
				+ l[c] * (dl_dx[a][b] + dl_dx[b][a])
			)
	end)
	-- conn^a_bc = g^ad conn_dbc
	return gU * connL
end

local metrics = table{
	{vacuum = Conn_Vacuum},
	{Newton = Conn_Newton},
	{Schwarzschild = Conn_Schwarzschild},
	{['Kerr-Schild'] = Conn_KerrSchild},
}

local metricNames = metrics:map(function(p) return (next(p)) end)

local currentMetric = ffi.new('int[1]', 2)
local Conn = select(2, next(metrics[currentMetric[0]+1]))

local function deriv(pt)
	local vm = matrix{pt.vel:unpack()}
	return pt_t(
		pt.vel,
		-vec3d((Conn(pt.pos) * vm * vm):unpack())
	)
end

_G.angleDistThreshold = 1
_G.angleThreshold = 1
_G.distThreshold = 1
_G.dt = .01
function App:simulate()
	if math.floor(self.t - dt) ~= math.floor(self.t) then
		rings:insert(pts:map(function(pt)
			return pt_t(pt.pos)
		end))
	end

	-- update by velocity
	self.t = self.t + dt
	for i=1,#pts do
		local pt = pts[i]
	
		-- [[ Euler
		local delta = deriv(pt) * dt
		--]]
		--[[ Runge-Kutta 4
		local k1 = deriv(pt) * dt 
		local k2 = deriv(pt + k1 * .5) * dt
		local k3 = deriv(pt + k2 * .5) * dt
		local k4 = deriv(pt + k3) * dt
		local delta = (k1 + k2 * 2 + k3 * 2 + k4) * (1/6) 
		--]]
		
		delta = delta * (dt / delta.pos.x)
		pts[i] = pts[i] + delta
		--pt.vel = pt.vel / math.sqrt(pt.vel[1]*pt.vel[1] - pt.vel[2] * pt.vel[2] - pt.vel[3] * pt.vel[3])
		
		--[[ re-normalize
		pt.vel.x = 1
		local beta = math.sqrt(pt.vel.x * pt.vel.x + pt.vel.y * pt.vel.y)
		pt.vel.y = pt.vel.y / beta
		pt.vel.z = pt.vel.z / beta
		--]]
	end

	-- measure angle range
	minAngle = math.huge
	maxAngle = -math.huge
	minDist = math.huge
	maxDist = -math.huge

	local i = 1
	local divs
	while i <= #pts do
		local ia = i
		local ib = ia % #pts + 1
		local ic = ib % #pts + 1
		local pa = pts[ia]
		local pb = pts[ib]
		local pc = pts[ic]

		local d1 = pb.pos - pa.pos
		local d1len = math.sqrt(d1.y * d1.y + d1.z * d1.z)
		d1 = d1 / d1len
		
		local d2 = pc.pos - pb.pos
		local d2len = math.sqrt(d2.y * d2.y + d2.z * d2.z)
		d2 = d2 / d2len
	
		local angle = math.acos(math.clamp(d1.y * d2.y + d1.z * d2.z, -1, 1))
		minAngle = math.min(minAngle, angle)
		maxAngle = math.max(maxAngle, angle)
		minDist = math.min(minDist, d1len)
		maxDist = math.max(maxDist, d1len)

		-- then subdivide

		if d1len * (angle + 1) > angleDistThreshold then
		--if d1len > distThreshold or angle > angleThreshold then
			if not divs then divs = {} end
			divs[ia] = true
		end
		if d2len * (angle + 1) > angleDistThreshold then
		--if d2len > distThreshold or angle > angleThreshold then
			if not divs then divs = {} end
			divs[ib] = true
		end
		
		i = i + 1
	end
	if divs then
		for i=#pts,1,-1 do
			if divs[i] then
				local pa = pts[i]
				local pb = pts[i%#pts+1]
				local avel = math.sqrt(pa.vel.y * pa.vel.y + pa.vel.z * pa.vel.z)
				local bvel = math.sqrt(pb.vel.y * pb.vel.y + pb.vel.z * pb.vel.z)
				
				local vel = (pa.vel + pb.vel) * .5
				
				vel.x = 1
				local beta = math.sqrt(vel.y * vel.y + vel.z * vel.z)
				vel.y = vel.y / beta
				vel.z = vel.z / beta
			
				-- TODO spherical average.  trace back from pa,pb pos along vel to find intersect in xy, then do spherical averaging
				local np = pt_t(
					(pa.pos + pb.pos) * .5,
					vel
				)
				pts:insert(i+1, np)
			end
		end
	end
end

function App:event(event, eventPtr)
	App.super.event(self, event, eventPtr)
	local canHandleMouse = not ig.igGetIO()[0].WantCaptureMouse
	local canHandleKeyboard = not ig.igGetIO()[0].WantCaptureKeyboard
	if event.type == sdl.SDL_KEYDOWN or event.type == sdl.SDL_KEYUP then
		if canHandleKeyboard and event.type == sdl.SDL_KEYDOWN then
			if event.key.keysym.sym == sdl.SDLK_SPACE then
				self.updateMethod = not self.updateMethod
			elseif event.key.keysym.sym == ('u'):byte() then
				self.updateMethod = 'step'
			elseif event.key.keysym.sym == ('r'):byte() then
				print'resetting...'
				self:reset()
				self.updateMethod = nil
			end
		end
	end
end

function App:update(...)
	gl.glClear(bit.bor(gl.GL_COLOR_BUFFER_BIT, gl.GL_DEPTH_BUFFER_BIT))

	gl.glMatrixMode(gl.GL_MODELVIEW)
	gl.glScalef(1,1,1)

	gl.glEnable(gl.GL_BLEND)
	gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE)

	gl.glColor3f(0,1,0)
	gl.glBegin(gl.GL_LINES)
	gl.glVertex3f(body.centerX, body.centerY, 0)
	gl.glVertex3f(body.centerX, body.centerY, self.t)
	gl.glEnd()
	for ir=0,1 do
		local r = ir == 0 and body.R or body.rs
		for it=0,1 do
			local t = it * self.t
			gl.glBegin(gl.GL_LINE_LOOP)
			local n = 20
			for i=1,n do
				local th = (i-.5)/n * 2 * math.pi
				gl.glVertex3d(r * math.cos(th) + body.centerX, r * math.sin(th) + body.centerY, t)
			end
			gl.glEnd()
		end
	end

	gl.glColor3f(1,1,1)
	gl.glBegin(gl.GL_LINE_LOOP)
	for _,pt in ipairs(pts) do
		gl.glVertex3d(pt.pos.y, pt.pos.z, pt.pos.x)
	end
	gl.glEnd()

	gl.glPointSize(3)
	gl.glBegin(gl.GL_POINTS)
	for _,pt in ipairs(pts) do
		gl.glVertex3d(pt.pos.y, pt.pos.z, pt.pos.x)
	end
	gl.glEnd()

	for i=1,#rings do
		local rowa = rings[i]
		local rowb = rings[i+1] or pts
		local na = #rowa
		local nb = #rowb
		local m = math.max(na,nb)
		gl.glColor3f(1,1,0)
		do local mode = gl.GL_FILL
		--for _,mode in ipairs{gl.GL_FILL, gl.GL_LINE} do
			if mode == gl.GL_FILL then
				gl.glColor4f(1,1,0,.1)
			else
				gl.glColor3f(1,1,1)
			end
			--gl.glPolygonMode(gl.GL_FRONT_AND_BACK, mode)
			gl.glBegin(gl.GL_TRIANGLE_STRIP)
			for j=0,m do
				local ia = math.floor((j+.5)/m*na)%na + 1
				local ib = math.floor((j+.5)/m*nb)%nb + 1
				gl.glVertex3d(rowa[ia].pos.y, rowa[ia].pos.z, rowa[ia].pos.x)
				gl.glVertex3d(rowb[ib].pos.y, rowb[ib].pos.z, rowb[ib].pos.x)
			end
			gl.glEnd()
		end
		--gl.glPolygonMode(gl.GL_FRONT_AND_BACK, gl.GL_FILL)
		gl.glColor3f(1,1,1)
		gl.glBegin(gl.GL_LINE_LOOP)
		for ia=1,na do
			gl.glVertex3d(rowa[ia].pos.y, rowa[ia].pos.z, rowa[ia].pos.x)
		end
		gl.glEnd()
	end
	gl.glColor3f(1,1,1)
	gl.glBegin(gl.GL_LINE_LOOP)
	for i,pt in ipairs(pts) do
		gl.glVertex3d(pt.pos.y, pt.pos.z, pt.pos.x)
	end
	gl.glEnd()

	if self.updateMethod then 
		if self.updateMethod == 'step' then self.updateMethod = nil end
		self:simulate() 
	end

	App.super.update(self, ...)
end

local buffer = ffi.new('char[256]')
local function inputFloat(name, t, k)
	local s = tostring(t[k])
	ffi.copy(buffer, s, #s)
	buffer[#s] = 0
	if ig.igInputText(name, buffer, ffi.sizeof(buffer)) then
		local v = tonumber(ffi.string(buffer))
		if v then
			t[k] = v
		end
	end
end

function App:updateGUI()
	ig.igText('vertex count '..#pts)
	ig.igText('time '..self.t)
	ig.igText('min angle '..minAngle)
	ig.igText('max angle '..maxAngle)
	ig.igText('min dist '..minDist)
	ig.igText('max dist '..maxDist)
	inputFloat('dt', _G, 'dt')
	inputFloat('angle threshold', _G, 'angleDistThreshold')
	--inputFloat('angle threshold', _G, 'angleThreshold')
	--inputFloat('dist threshold', _G, 'distThreshold')
	inputFloat('body center x', body, 'centerX')
	inputFloat('body center y', body, 'centerY')
	inputFloat('body Schwarzschild radius', body, 'R')
	inputFloat('body surface radius', body, 'rs')
	if ig.igCombo('metric', currentMetric, metricNames) then
		Conn = select(2, next(metrics[currentMetric[0]+1]))
	end
end

App():run()
