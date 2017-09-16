#include <DDImage/Iop.h>
#include <DDImage/Knobs.h>
#include <DDImage/Row.h>
#include <DDImage/Tile.h>
#include <DDImage/Thread.h>
#include <memory>
#include "nlmeans_lib.h"

using namespace DD::Image;

static const char* const CLASS = "Nlmeans";
static const char* const HELP = "Denoises an image using the Buades/Coll/Morel non-local means algorithm";

class Nlmeans : public Iop {
	private:
		float 	sigma, multiplier;
		float	bloc;
		float	*wholeframe_in, *wholeframe_out;
		int	rx, ry, rr, rt;
		Lock	bufferlock;
	public:
				Nlmeans(Node*);
				~Nlmeans();
		
		const char* 	Class() const { return CLASS; }
		const char* 	node_help() const { return HELP; }
		static const 	Description d;

		void 		_validate(bool);
		void 		_request(int, int, int, int, ChannelMask, int);
		void 		engine(int, int, int, ChannelMask, Row&);
		
		virtual void 	knobs(Knob_Callback);
};

static Iop		*Nlmeans_create(Node *n) { return new Nlmeans(n); }
const Iop::Description 	Nlmeans::d(CLASS, "Filter/Nlmeans", Nlmeans_create);

Nlmeans::Nlmeans(Node *n) : Iop(n) {
	sigma = 0.01;
	multiplier = 1.5;
	bloc = 5.0;
	wholeframe_in = NULL;
	wholeframe_out = NULL;
}

Nlmeans::~Nlmeans() {
	if(wholeframe_in != NULL) {
		free(wholeframe_in);
		wholeframe_in = NULL;
		free(wholeframe_out);
		wholeframe_out = NULL;
	}
}

void Nlmeans::knobs(Knob_Callback kcb) {
	Float_knob(kcb, &sigma, IRange(0,0.1), "sigma");
	Tooltip(kcb, "The standard deviation of the noise to be removed: more = harder");
	//Float_knob(kcb, &multiplier, IRange(0,3), "multiplier");
	Float_knob(kcb, &bloc, IRange(1,15), "bloc");
	Tooltip(kcb, "The size of the search area: more = slower");
}

void Nlmeans::_validate(bool for_real) {
	//std::cout << "Nlmeans: in _validate()" << endl;
		
	copy_info();
		
	ChannelSet cs;
	cs.clear();
	cs += Chan_Red;
	cs += Chan_Green;
	cs += Chan_Blue;
	set_out_channels(cs);
}

void Nlmeans::_request(int x, int y, int r, int t, ChannelMask channels, int count) {
	//std::cout << "Nlmeans: in _request()" << endl;
		
	if(wholeframe_in != NULL) {
		free(wholeframe_in);
		wholeframe_in = NULL;
		free(wholeframe_out);
		wholeframe_out = NULL;
	}
	
	rx = x;
	ry = y;
	rr = r;
	rt = t;
	
	input(0)->new_request_pass();
	input(0)->request(x, y, r, t, channels, count);
	
	//std::cout << "Nlmeans: got _request for " << x << " " << y << " " << r << " " << t << endl;
}

void Nlmeans::engine(int y, int x, int r, ChannelMask channels, Row &out) {
	// std::cout << "Nlmeans: in engine()" << endl;
	
	Tile	t(input0(), rx, ry, rr, rt, channels);
 	if(aborted()) return;
	t.load_range(t.y(), t.t());
	
	bufferlock.lock();
	if(wholeframe_in == NULL) {	
		// Our private buffer is empty, fill it from Nuke's tile
		wholeframe_in = (float *)std::malloc(t.w() * t.h() * 3 * sizeof(float));
		wholeframe_out = (float *)std::malloc(t.w() * t.h() * 3 * sizeof(float));
		if(wholeframe_in == NULL || wholeframe_out == NULL) {
			std::cout << "Nlmeans: out of memory!" << endl;
			return;
		} else {
			//std::cout << "Nlmeans: allocated buffers for " << t.w() << " x " << t.h() << " x " << 3 << endl;
		}
		for(int i = 0; i < t.h(); i++) {
			std::memcpy(&wholeframe_in[i * t.w()], &t[Chan_Red][t.y() + i][t.x()], t.w() * sizeof(float));
			std::memcpy(&wholeframe_in[(t.h() * t.w()) + (i * t.w())], &t[Chan_Green][t.y() + i][t.x()], t.w() * sizeof(float));
			std::memcpy(&wholeframe_in[(t.h() * t.w() * 2) + (i * t.w())], &t[Chan_Blue][t.y() + i][t.x()], t.w() * sizeof(float));
		}
		//std::cout << "Nlmeans: filled buffer" << endl;
		
		
		// Woohoo! Denoise time.
        	int win = 1;
        	int averageoption = 1;
        	float *kernel = NULL;
        	int ksize = -1;
		//std::cout << "Nlmeans: denoising..." << flush;
        	nlmeans(win, bloc, averageoption, kernel, ksize, multiplier, sigma, sigma, sigma,
				wholeframe_in, wholeframe_in + t.w() * t.h(), wholeframe_in + 2 * t.w() * t.h(),
                		wholeframe_out, wholeframe_out + t.w() * t.h(), wholeframe_out + 2 * t.w() * t.h(), t.w(), t.h(), this);
		//std::cout << " done." << endl;
	}		
	bufferlock.unlock();
	
	// Now, we've got a whole buffer ready. Feed Nuke the pixels it wants.
	foreach(channel, channels) {
		float *ptr = out.writable(channel);
		// std::cout << "Nlmeans: filling channel " << channel << ", y=" << y << ", x from " << x << " to " << r << endl;
		switch(channel) {
			case Chan_Red:
				std::memcpy(&ptr[x], &wholeframe_out[((y - t.y()) * t.w()) - t.x() + x], (r - x) * sizeof(float));
				break;
			case Chan_Green:
				std::memcpy(&ptr[x], &wholeframe_out[(t.h() * t.w()) + ((y - t.y()) * t.w()) - t.x() + x], (r - x) * sizeof(float));
				break;
			case Chan_Blue:
				std::memcpy(&ptr[x], &wholeframe_out[(t.h() * t.w() * 2) + ((y - t.y()) * t.w()) - t.x() + x], (r - x) * sizeof(float));
				break;
			default:
				std::memcpy(&ptr[x], &t[channel][y][x], (r - x) * sizeof(float));
		}
	}
		
}













